"""健康计算器 (/health) —— zhuyidao.net 姊妹应用的字节级复刻.

输入校验顺序、公式、阈值、状态文案与提示语全部来自线上实测
（bench/data/health-probes.jsonl，规则推导见 docs/health-notes.md）。
入口 ``validate(data)`` 返回 ``(fields, err)``、``calculate(fields)`` 返回
``results_dict``；server 层负责
把错误包成 ``{"error": ..., "ok": false}``、把结果包成
``{"ok": true, "record_id": n, "results": ...}``。
"""

import contextlib
import json
import math
import os
import re
import threading
from pathlib import Path

# ---- 校验 -----------------------------------------------------------------

# 昵称字符集：字母/数字/下划线/短横线/中文（\w 的 Unicode 类覆盖了 café、
# 日文、希腊文等情况——由探针 s2:nick-jp / s2:nick-greek 钉死）
NICK_RE = re.compile(r"^[\w一-鿿-]+$")

PAL_CHOICES = (1.55, 1.85, 2.2)

EXERCISE_TYPES = {
    "walking_slow": ("慢速散步", 2.8),
    "walking_moderate": ("普通速度走路", 3.8),
    "walking_brisk": ("快走", 4.8),
    "jogging": ("慢跑", 7.5),
    "running_8kph": ("跑步（约8公里/小时）", 8.5),
    "cycling_easy": ("轻松骑车", 4.3),
    "cycling_moderate": ("中等速度骑车", 7.0),
    "swimming_leisure": ("休闲游泳", 6.0),
    "strength_training": ("一般力量训练", 3.5),
    "yoga": ("普通瑜伽", 2.3),
    "badminton": ("休闲羽毛球", 5.5),
    "table_tennis": ("乒乓球", 4.0),
    "stair_climbing": ("一般速度爬楼梯", 6.8),
    "jump_rope": ("跳绳", 11.0),
}

BP_CONTEXTS = {"clinic": "诊室测量", "home": "家庭自测"}

# CUN-BAE / Deurenberg 的“结果超出公式的合理解释范围”界值（探针钉死：
# 两公式共用 [0,70]——bf -0.1/70.6 触发、0.3/69.5 不触发；deur -0.1/70.3
# 触发、0.5/68.7 不触发。超龄状态不豁免，cascade 只看数值旗标）
BF_FLAG_LO, BF_FLAG_HI = 0.0, 70.0

# 身高提示触发界值：h=99 触发、h=100 不触发 -> h<100
HEIGHT_TIP_CM = 100.0


def missing(v):
    """None 或空字符串视为未填（前端对留空字段发 null；'' 实测同 absent）。"""
    return v is None or v == ""


def num(v, label):
    """站端数值解析，三类错误（s2/s3 实测）：bool -> '格式不正确'；
    float() 失败 -> '必须是数字'；NaN/±Inf -> '必须是有限数字'。"""
    if isinstance(v, bool):
        return None, f"{label}格式不正确。"
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None, f"{label}必须是数字。"
    if not math.isfinite(x):
        return None, f"{label}必须是有限数字。"
    return x, None


def num_field(data, key, label, lo, hi, required=False, integer=False):
    """数值字段校验：required 缺席 -> 请填写{label}；选填缺席 -> (None, None)。
    integer 在范围之后追加整数检查并返回 int（范围先于整数，探针钉死）。"""
    v = data.get(key)
    if missing(v):
        return (None, f"请填写{label}。") if required else (None, None)
    x, err = num(v, label)
    if err:
        return None, err
    if not lo <= x <= hi:
        return None, f"{label}应在 {lo}～{hi} 之间。"
    if integer:
        if not x.is_integer():
            return None, f"{label}必须是整数。"
        x = int(x)
    return x, None


def validate(data):
    """按站端顺序校验请求体；返回 (cleaned_fields, error_message)。

    顺序（成对探针钉死）：蜜罐 -> 昵称 -> 性别 -> 年龄 -> 身高 -> 体重 ->
    活动水平 -> 腰围 -> 臀围 -> 静息心率 -> 运动(成对) -> 血压(成组) -> 睡眠。
    """
    f = {}
    # 蜜罐字段：任何真值都拒绝
    if data.get("website"):
        return None, "提交未通过校验。"

    # 昵称：str() 强转 + trim；空 -> 请填写；超长 -> 长度错（先于字符集，
    # s2:nick-25-bad 实测）；非法字符 -> 字符集错
    nick_raw = data.get("nickname")
    nick = str(nick_raw).strip() if nick_raw is not None else ""
    if not nick:
        return None, "请填写昵称。"
    if len(nick) > 24:
        return None, "昵称不能超过24个字符。"
    if not NICK_RE.match(nick):
        return None, "昵称只能包含中文、字母、数字、下划线或短横线。"
    f["nickname"] = nick

    # 性别：大小写归一（"MALE" -> male 实测通过）
    sex_raw = data.get("sex")
    sex = str(sex_raw).strip().lower() if sex_raw is not None else ""
    if sex not in ("male", "female"):
        return None, "请选择性别，以便使用对应公式。"
    f["sex"] = sex

    # 年龄：数字 -> 范围(18~120) -> 整数（17.9 实测报范围错，范围先于整数检查）
    x, err = num_field(data, "age", "年龄", 18, 120, required=True, integer=True)
    if err:
        return None, err
    f["age"] = x

    x, err = num_field(data, "height_cm", "身高", 1, 250, required=True)
    if err:
        return None, err
    f["height_cm"] = x

    x, err = num_field(data, "weight_kg", "体重", 1, 500, required=True)
    if err:
        return None, err
    f["weight_kg"] = x

    # 活动水平：数字 -> [1.4,2.4] -> 三档枚举
    x, err = num_field(data, "pal", "活动水平", 1.4, 2.4, required=True)
    if err:
        return None, err
    if x not in PAL_CHOICES:
        return None, "请选择页面提供的活动水平。"
    f["pal"] = x

    # 选填数值，顺序：腰围 -> 臀围 -> 静息心率 -> 运动 -> 血压 -> 睡眠
    x, err = num_field(data, "waist_cm", "腰围", 1, 300)
    if err:
        return None, err
    f["waist_cm"] = x

    x, err = num_field(data, "hip_cm", "臀围", 1, 300)
    if err:
        return None, err
    f["hip_cm"] = x

    x, err = num_field(data, "resting_hr", "静息心率", 20, 250, integer=True)
    if err:
        return None, err
    f["resting_hr"] = x

    # 运动：类型在场先查枚举（s2:extype-bogus-solo 实测），再查成对，再查时长
    etype, emin = data.get("exercise_type"), data.get("exercise_minutes")
    if not missing(etype) and str(etype) not in EXERCISE_TYPES:
        return None, "请选择页面提供的运动类型。"
    if missing(etype) != missing(emin):
        return None, "运动类型和运动时长需要一起填写，或都留空。"
    if not missing(etype):
        f["exercise_type"] = str(etype)
        x, err = num_field(data, "exercise_minutes", "运动时长", 1, 1440, integer=True)
        if err:
            return None, err
        f["exercise_minutes"] = x
    else:
        f["exercise_type"] = None
        f["exercise_minutes"] = None

    # 血压：场景在场先查枚举（s2:bpctx-bogus-solo 实测），再查成组，
    # 再收缩压 -> 舒张压 -> 收缩压须高于舒张压（s3:bp-40-95/100-100 钉死顺序）
    bctx, bsys, bdia = (data.get("bp_context"), data.get("systolic_bp"), data.get("diastolic_bp"))
    if not missing(bctx) and str(bctx) not in BP_CONTEXTS:
        return None, "请选择页面提供的血压测量场景。"
    present = [not missing(x) for x in (bctx, bsys, bdia)]
    if any(present) and not all(present):
        return None, "血压测量场景、收缩压和舒张压需要一起填写，或都留空。"
    if all(present):
        s, err = num_field(data, "systolic_bp", "收缩压", 50, 300, integer=True)
        if err:
            return None, err
        d, err = num_field(data, "diastolic_bp", "舒张压", 30, 200, integer=True)
        if err:
            return None, err
        if s <= d:
            return None, "收缩压通常应高于舒张压，请检查输入。"
        f["bp_context"], f["systolic_bp"], f["diastolic_bp"] = str(bctx), s, d
    else:
        f["bp_context"] = f["systolic_bp"] = f["diastolic_bp"] = None

    x, err = num_field(data, "sleep_hours", "平均每晚睡眠时长", 0, 24)
    if err:
        return None, err
    f["sleep_hours"] = x

    return f, None


# ---- 结果组装 -------------------------------------------------------------

NHC_WEIGHT_URL = "https://www.nhc.gov.cn/wjw/ylyjs/202412/b3d40e0141834897808ce6c9dce76a60.shtml"
NHC_WAIST_URL = "https://www.nhc.gov.cn/xcs/c100122/202503/27b26d397f4941ff8ff3846ef40026ee.shtml"
NHC_BP_URL = "https://www.nhc.gov.cn/cms-search/downFiles/63f752a17cfd4b4781f744477561866f.pdf"
NHC_SLEEP_URL = (
    "https://www.nhc.gov.cn/guihuaxxs/c100133/202503/70d5836afe804a858b899ee951a24a13.shtml"
)
AHA_HR_URL = (
    "https://www.heart.org/en/healthy-living/"
    "exercise-and-physical-activity/fitness-basics/"
    "target-heart-rates"
)
CDC_KARVONEN_URL = "https://stacks.cdc.gov/view/cdc/86801"
CUNBAE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC3263863/"
DEURENBERG_URL = "https://pubmed.ncbi.nlm.nih.gov/2043597/"
WATSON_URL = "https://pubmed.ncbi.nlm.nih.gov/6986753/"


def item(
    key,
    category,
    name,
    display,
    numeric,
    unit,
    status,
    direction,
    reference,
    sort_order,
    url,
    metadata=None,
    version="1.0",
):
    """组装一个指标卡片；direction 非空时同步写进 metadata（站端行为）。"""
    meta = dict(metadata) if metadata else {}
    if direction is not None:
        meta["direction"] = direction
    return {
        "key": key,
        "category_key": category,
        "sort_order": sort_order,
        "display_name": name,
        "display_value": display,
        "numeric_value": numeric,
        "unit": unit,
        "status": status,
        "direction": direction,
        "reference": reference,
        "source_url": url,
        "formula_version": version,
        "metadata": meta,
    }


def num_item(
    key,
    category,
    name,
    value,
    digits,
    unit,
    status,
    direction,
    reference,
    sort_order,
    url,
    version="1.0",
):
    """数值型卡片：numeric/display 均为 round(value, digits)。"""
    rounded = round(value, digits)
    return item(
        key,
        category,
        name,
        str(rounded),
        rounded,
        unit,
        status,
        direction,
        reference,
        sort_order,
        url,
        version=version,
    )


def range_item(key, category, name, low, high, digits, unit, reference, sort_order, url):
    """区间型卡片：display 'low～high'，numeric 为 null，metadata 带边界。"""
    lo, hi = round(low, digits), round(high, digits)
    return item(
        key,
        category,
        name,
        f"{lo}～{hi}",
        None,
        unit,
        None,
        None,
        reference,
        sort_order,
        url,
        metadata={"low": lo, "high": hi},
    )


# ---- 提示语（中英成对；触发顺序见 docs/health-notes.md）-------------------


def build_tips(
    fields,
    bmi_raw,
    bmi_status,
    wadj,
    waist_status,
    bf_flag,
    bf_out_of_age,
    watson_flag,
    rhr_status,
    bp_direction,
    sleep_dir,
):
    """按站端固定顺序生成 tips / tips_en。"""
    tips, tips_en = [], []
    if fields["height_cm"] < HEIGHT_TIP_CM:
        tips.append(
            "身高明显低于成人常见范围；虽然通过了1～250厘米的硬校验，但这些成人公式很可能不适用。"
        )
        tips_en.append(
            "The height is well below the range commonly seen in "
            "adults. Although it passes the hard 1–250 cm "
            "validation, these adult formulas are unlikely to "
            "apply."
        )
    if wadj is not None:
        amount, gain = wadj
        if gain:
            tips.append(
                f"按BMI正常范围反算，至少还需增重{amount}公斤可达到"
                "BMI 18.5；这只是数学反算，不代表个人治疗目标。"
            )
            tips_en.append(
                f"Based on the BMI reference range, gaining at "
                f"least {amount} kg would reach a BMI of 18.5. "
                "This is only a mathematical conversion, not an "
                "individual treatment target."
            )
        else:
            tips.append(
                f"按BMI正常范围反算，至少还需减重{amount}公斤可回到"
                "BMI 23.9；这只是数学反算，不代表个人治疗目标。"
            )
            tips_en.append(
                f"Based on the BMI reference range, losing at "
                f"least {amount} kg would return BMI to 23.9. "
                "This is only a mathematical conversion, not an "
                "individual treatment target."
            )
    if bmi_raw < 18.5:
        tips.append("BMI低于18.5；若近期体重持续下降、食欲差或容易疲劳，建议咨询专业人员。")
        tips_en.append(
            "BMI is below 18.5. Consider professional advice if "
            "weight has continued to fall or if poor appetite or "
            "fatigue is present."
        )
    elif bmi_raw >= 24:
        tips.append("BMI达到超重或肥胖参考范围；建议优先采用可持续的饮食与运动调整，避免快速节食。")
        tips_en.append(
            "BMI is in the overweight or obesity reference range. "
            "Prefer sustainable diet and activity changes and "
            "avoid rapid dieting."
        )
    else:
        tips.append("BMI处于中国成人正常参考范围，可继续结合腰围、运动习惯和身体成分观察趋势。")
        tips_en.append(
            "BMI is within the Chinese reference range for adults. "
            "Continue to consider waist circumference, activity "
            "habits, body composition, and long-term trends."
        )
    if waist_status in ("中心型肥胖前期", "达到中心型肥胖界值"):
        tips.append("腰围提示腹部脂肪风险升高；可定期在同一位置复测，并关注血压、血糖和血脂。")
        tips_en.append(
            "Waist circumference suggests increased abdominal-fat "
            "risk. Recheck periodically at the same site and also "
            "monitor blood pressure, glucose, and lipids."
        )
    # 超龄提示与不宜提示并存（s3:cb-a120-* 实测：a120+不宜 时两条都出）
    if bf_flag:
        tips.append("CUN-BAE给出了不宜解释的结果，身体成分与Cunningham公式已停止计算。")
        tips_en.append(
            "CUN-BAE produced a result that should not be "
            "interpreted, so body-composition and Cunningham "
            "calculations were stopped."
        )
    else:
        tips.append(
            "CUN-BAE作为主要体脂估算，Deurenberg作为对照；两者都是"
            "统计公式，不能替代身体成分测量或医学诊断。"
        )
        tips_en.append(
            "CUN-BAE is the primary body-fat estimate and "
            "Deurenberg is shown for comparison. Both are "
            "statistical formulas and cannot replace "
            "body-composition measurement or medical diagnosis."
        )
    if bf_out_of_age:
        tips.append("CUN-BAE原始研究对象为18～80岁；当前年龄超出其验证范围，结果仅供查看。")
        tips_en.append(
            "The original CUN-BAE study covered ages 18–80. "
            "The entered age is outside its validation range, "
            "so the result is for reference only."
        )
    if watson_flag:
        tips.append("Watson体水分公式给出了不宜解释的结果，因此体水分及其百分比未展示。")
        tips_en.append(
            "The Watson formula produced an uninterpretable "
            "result, so total body water and its percentage are "
            "not shown."
        )
    if rhr_status is not None and rhr_status != "一般成人常见范围":
        tips.append(
            "静息心率受运动水平、情绪、药物和疾病影响；若多次异常或伴有不适，应咨询医务人员。"
        )
        tips_en.append(
            "Resting heart rate is affected by fitness, emotion, "
            "medication, and illness. Seek medical advice if "
            "readings are repeatedly unusual or accompanied by "
            "symptoms."
        )
    if bp_direction == "up":
        tips.append(
            "血压读数高于理想范围；单次测量不能诊断高血压，建议按"
            "规范在不同时间复测，持续异常时咨询医务人员。"
        )
        tips_en.append(
            "The blood-pressure reading is above the ideal range. "
            "One reading cannot diagnose hypertension; repeat "
            "measurements correctly at different times and seek "
            "medical advice if they remain abnormal."
        )
    elif bp_direction == "down":
        tips.append(
            "血压读数处于偏低范围；若反复出现或伴有头晕、乏力、视物模糊或晕厥，建议咨询医务人员。"
        )
        tips_en.append(
            "The blood-pressure reading is in a low range. Seek "
            "medical advice if this recurs or is accompanied by "
            "dizziness, weakness, blurred vision, or fainting."
        )
    if sleep_dir == "down":
        tips.append(
            "平均睡眠时长低于参考范围；除时长外，也应结合醒后精神"
            "状态、夜间觉醒和作息规律判断睡眠质量。"
        )
        tips_en.append(
            "Average sleep duration is below the reference range. "
            "Also consider how you feel after waking, night-time "
            "awakenings, and sleep regularity."
        )
    elif sleep_dir == "up":
        tips.append(
            "平均睡眠时长高于参考范围；若长期睡得较多仍感到疲倦，建议关注睡眠质量并咨询专业人员。"
        )
        tips_en.append(
            "Average sleep duration is above the reference range. "
            "If long sleep persists while fatigue remains, "
            "consider sleep quality and seek professional advice."
        )
    return tips, tips_en


def body_fat_status(value, male, flagged):
    """体脂率 (status, direction) 分档：flagged -> 不宜解释；否则按性别界值分档。"""
    if flagged:
        return "结果超出公式的合理解释范围", None
    warn = 20.0 if male else 25.0
    alarm = 25.0 if male else 30.0
    if value < warn:
        return "未达到高体脂报警界值", None
    if value < alarm:
        return "处于体脂偏高参考区间", "up"
    return "达到高体脂报警界值", "up"


def calculate(f):
    """对校验通过的字段求值，返回 results 字典（站端 wire 结构）。"""
    sex = f["sex"]
    male = sex == "male"
    age, h, w = f["age"], f["height_cm"], f["weight_kg"]
    h_m = h / 100
    bmi_raw = w / h_m**2
    bmi = round(bmi_raw, 1)
    waist, hip = f["waist_cm"], f["hip_cm"]
    rhr = f["resting_hr"]

    # --- 体重 / BMI ---
    if bmi_raw < 18.5:
        bmi_status, bmi_dir = "体重过低", "down"
    elif bmi_raw < 24:
        bmi_status, bmi_dir = "正常范围", None
    elif bmi_raw < 28:
        bmi_status, bmi_dir = "超重范围", "up"
    else:
        bmi_status, bmi_dir = "肥胖范围", "up"
    tw_low = round(18.5 * h_m**2, 1)
    tw_high = round(23.9 * h_m**2, 1)
    wadj = None
    if bmi_raw < 18.5:
        wadj = (round(18.5 * h_m**2 - w, 1), True)
    elif bmi_raw >= 24:
        wadj = (round(w - 23.9 * h_m**2, 1), False)

    # --- 体脂：CUN-BAE（主）与 Deurenberg（对照），同一套界值 ---
    bf = (
        -44.988
        + 0.503 * age
        + 3.172 * bmi_raw
        - 0.026 * bmi_raw**2
        - 0.02 * bmi_raw * age
        + 0.00021 * bmi_raw**2 * age
    )
    if not male:
        bf += 10.689 + 0.181 * bmi_raw - 0.005 * bmi_raw**2
    deur = 1.2 * bmi_raw + 0.23 * age - (10.8 if male else 0) - 5.4
    bf_flag = not BF_FLAG_LO <= bf <= BF_FLAG_HI
    deur_flag = not BF_FLAG_LO <= deur <= BF_FLAG_HI
    bf_out_of_age = age > 80
    if bf_out_of_age:
        bf_status, bf_dir = "超出CUN-BAE原始18～80岁验证范围", None
    else:
        bf_status, bf_dir = body_fat_status(bf, male, bf_flag)
    deur_status, deur_dir = body_fat_status(deur, male, deur_flag)

    # CUN-BAE 不宜 -> 派生体成分与 Cunningham 停止（超龄不停）
    fat_mass = round(bf * w / 100, 1) if not bf_flag else None
    ffm = round(w - bf * w / 100, 1) if not bf_flag else None
    janma = 9270 * w / (6680 + 216 * bmi_raw) if male else 9270 * w / (8780 + 244 * bmi_raw)
    bsa = math.sqrt(h * w / 3600)

    # --- 腰围 / 腰高比 / 腰臀比 / 派生体型指数 ---
    if waist is not None:
        if male:
            waist_status = (
                "未达到中心型肥胖前期界值"
                if waist < 85
                else "中心型肥胖前期"
                if waist < 90
                else "达到中心型肥胖界值"
            )
        else:
            waist_status = (
                "未达到中心型肥胖前期界值"
                if waist < 80
                else "中心型肥胖前期"
                if waist < 85
                else "达到中心型肥胖界值"
            )
        waist_dir = None if waist_status == "未达到中心型肥胖前期界值" else "up"
        whtr_raw = waist / h
        whtr = round(whtr_raw, 3)
        whtr_status = "低于0.5参考界值" if whtr_raw < 0.5 else "达到或超过0.5参考界值"
        whtr_dir = None if whtr_raw < 0.5 else "up"
        rfm = round(64 - 20 * h / waist + (0 if male else 12), 1)
        absi = round((waist / 100) / (bmi_raw ** (2 / 3) * math.sqrt(h_m)), 4)
        bri = round(
            364.2 - 365.5 * math.sqrt(1 - ((waist / 100 / (2 * math.pi)) ** 2) / (0.5 * h_m) ** 2),
            2,
        )
        conicity = round((waist / 100) / (0.109 * math.sqrt(w / h_m)), 3)
    else:
        waist_status = waist_dir = whtr = whtr_status = whtr_dir = None
        rfm = absi = bri = conicity = None
    bai = round(hip / h_m**1.5 - 18, 1) if hip is not None else None
    whr_limit = 0.9 if male else 0.85
    if waist is not None and hip is not None:
        whr_raw = waist / hip
        whr = round(whr_raw, 3)
        whr_status = "低于参考界值" if whr_raw < whr_limit else "达到或超过参考界值"
        whr_dir = None if whr_raw < whr_limit else "up"
    else:
        whr = whr_status = whr_dir = None

    # --- Watson 体水分 / Nadler 血容量 ---
    tbw = (
        2.447 - 0.09516 * age + 0.1074 * h + 0.3362 * w
        if male
        else -2.097 + 0.1069 * h + 0.2466 * w
    )
    watson_flag = tbw <= 0 or tbw > w
    watson_tbw = None if watson_flag else round(tbw, 1)
    water_pct = None if watson_flag else round(tbw / w * 100, 1)
    blood = (
        0.3669 * h_m**3 + 0.03219 * w + 0.6041 if male else 0.3561 * h_m**3 + 0.03308 * w + 0.1833
    )

    # --- 代谢 ---
    mifflin = 10 * w + 6.25 * h - 5 * age + (5 if male else -161)
    harris = (
        88.362 + 13.397 * w + 4.799 * h - 5.677 * age
        if male
        else 447.593 + 9.247 * w + 3.098 * h - 4.330 * age
    )
    cunningham = None if bf_flag else round(500 + 22 * (w - bf * w / 100), 0)
    tdee = round(mifflin * f["pal"], 0)

    # --- 心率 ---
    max_hr_raw = 208 - 0.7 * age
    max_hr = round(max_hr_raw, 0)
    mod_lo, mod_hi = round(max_hr_raw * 0.5, 0), round(max_hr_raw * 0.7, 0)
    vig_lo, vig_hi = round(max_hr_raw * 0.7, 0), round(max_hr_raw * 0.85, 0)
    if rhr is not None:
        if rhr < 60:
            rhr_status, rhr_dir = "低于一般成人常见范围", "down"
        elif rhr <= 100:
            rhr_status, rhr_dir = "一般成人常见范围", None
        else:
            rhr_status, rhr_dir = "高于一般成人常见范围", "up"
        hrr_raw = max_hr_raw - rhr
        hrr = round(hrr_raw, 0)
        karv_lo = round(rhr + 0.5 * hrr_raw, 0)
        karv_hi = round(rhr + 0.7 * hrr_raw, 0)
        vo2 = round(15.3 * max_hr_raw / rhr, 1)
    else:
        rhr_status = rhr_dir = hrr = karv_lo = karv_hi = vo2 = None

    # --- 营养 ---
    protein_lo, protein_hi = round(1.4 * w, 0), round(2.0 * w, 0)

    # --- 运动 ---
    if f["exercise_type"] is not None:
        label, met = EXERCISE_TYPES[f["exercise_type"]]
        kcal = round(met * 3.5 * w / 200 * f["exercise_minutes"], 0)
    else:
        label = met = kcal = None

    # --- 血压 / 睡眠 ---
    if f["bp_context"] is not None:
        ctx_label = BP_CONTEXTS[f["bp_context"]]
        sys_bp, dia_bp = f["systolic_bp"], f["diastolic_bp"]
        sys_status = (
            "低于90参考值" if sys_bp < 90 else "理想范围" if sys_bp < 120 else "高于理想范围"
        )
        sys_dir = "down" if sys_bp < 90 else "up" if sys_bp >= 120 else None
        dia_status = (
            "低于60参考值" if dia_bp < 60 else "理想范围" if dia_bp < 80 else "高于理想范围"
        )
        dia_dir = "down" if dia_bp < 60 else "up" if dia_bp >= 80 else None
        # 类别判断高侧优先（s3:bp-150-50 -> 1级 而非偏低；s3:bph-140-50 ->
        # 家庭筛查界值），偏低只在无任何高侧档位时兜底
        if f["bp_context"] == "home":
            if sys_bp >= 135 or dia_bp >= 85:
                bp_cat, bp_dir = "达到家庭血压高血压筛查界值", "up"
            elif sys_bp >= 120 or dia_bp >= 80:
                bp_cat, bp_dir = "正常高值范围", "up"
            elif sys_bp < 90 or dia_bp < 60:
                bp_cat, bp_dir = "血压偏低范围", "down"
            else:
                bp_cat, bp_dir = "正常血压范围", None
        else:
            if sys_bp >= 180 or dia_bp >= 110:
                bp_cat, bp_dir = "3级高血压范围", "up"
            elif sys_bp >= 160 or dia_bp >= 100:
                bp_cat, bp_dir = "2级高血压范围", "up"
            elif sys_bp >= 140 or dia_bp >= 90:
                bp_cat, bp_dir = "1级高血压范围", "up"
            elif sys_bp >= 120 or dia_bp >= 80:
                bp_cat, bp_dir = "正常高值范围", "up"
            elif sys_bp < 90 or dia_bp < 60:
                bp_cat, bp_dir = "血压偏低范围", "down"
            else:
                bp_cat, bp_dir = "正常血压范围", None
    else:
        ctx_label = sys_bp = dia_bp = sys_status = dia_status = None
        sys_dir = dia_dir = bp_cat = bp_dir = None

    if age >= 65:
        sleep_lo, sleep_hi, sleep_ref = 6.0, 7.0, "老年人一般6～7小时"
    else:
        sleep_lo, sleep_hi, sleep_ref = 7.0, 8.0, "成年人一般7～8小时"
    if f["sleep_hours"] is not None:
        s = f["sleep_hours"]
        if s < sleep_lo:
            sleep_status, sleep_dir = "低于睡眠时长参考范围", "down"
        elif s > sleep_hi:
            sleep_status, sleep_dir = "高于睡眠时长参考范围", "up"
        else:
            sleep_status, sleep_dir = "睡眠时长参考范围", None
    else:
        sleep_status = sleep_dir = None

    # --- 指标卡片（sort_order 即站端顺序）---
    items = [
        num_item(
            "bmi",
            "body",
            "BMI",
            bmi_raw,
            1,
            "kg/m²",
            bmi_status,
            bmi_dir,
            "中国成人18.5～23.9",
            10,
            "https://www.cdc.gov/bmi/about/index.html",
        ),
        range_item(
            "target_weight",
            "body",
            "BMI正常区间对应体重",
            tw_low,
            tw_high,
            1,
            "kg",
            "仅为BMI反算区间",
            20,
            NHC_WEIGHT_URL,
        ),
    ]
    if wadj is not None:
        amount, gain = wadj
        items.append(
            item(
                "weight_adjustment",
                "body",
                "达到BMI正常范围所需调整",
                str(amount),
                amount,
                "kg",
                "BMI体重过低，还需增重" if gain else "BMI达到超重或肥胖范围，还需减重",
                "down" if gain else "up",
                "增至BMI 18.5的数学反算值" if gain else "减至BMI 23.9的数学反算值",
                25,
                NHC_WEIGHT_URL,
                metadata={"adjustment_direction": "gain" if gain else "lose"},
            )
        )
    bf_ref = "男性20%～25%提示偏高，≥25%为报警值" if male else "女性25%～30%提示偏高，≥30%为报警值"
    items += [
        num_item(
            "body_fat_cun_bae",
            "body",
            "CUN-BAE估算体脂率",
            bf,
            1,
            "%",
            bf_status,
            bf_dir,
            bf_ref,
            30,
            CUNBAE_URL,
        ),
        num_item(
            "body_fat_deurenberg",
            "body",
            "Deurenberg对照体脂率",
            deur,
            1,
            "%",
            deur_status,
            deur_dir,
            "简化公式对照值",
            35,
            DEURENBERG_URL,
        ),
    ]
    if fat_mass is not None:
        items.append(
            num_item(
                "fat_mass",
                "body",
                "估算脂肪重量",
                fat_mass,
                1,
                "kg",
                None,
                None,
                "由估算体脂率推算",
                40,
                CUNBAE_URL,
                version="2.0",
            )
        )
    if ffm is not None:
        items.append(
            num_item(
                "ffm",
                "body",
                "估算去脂体重",
                ffm,
                1,
                "kg",
                None,
                None,
                "由估算体脂率推算",
                50,
                CUNBAE_URL,
                version="2.0",
            )
        )
    items += [
        num_item(
            "janma_lbm",
            "body",
            "Janmahasatian瘦体重",
            janma,
            1,
            "kg",
            None,
            None,
            "公式估算",
            60,
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3833312/",
        ),
        num_item(
            "bsa_mosteller",
            "body",
            "体表面积 BSA",
            bsa,
            2,
            "m²",
            None,
            None,
            "没有统一健康范围",
            70,
            "https://pubmed.ncbi.nlm.nih.gov/3657876/",
        ),
    ]
    if waist is not None:
        items += [
            num_item(
                "waist_assessment",
                "body",
                "腰围判断",
                waist,
                1,
                "cm",
                waist_status,
                waist_dir,
                "男性前期85、中心型90厘米" if male else "女性前期80、中心型85厘米",
                80,
                NHC_WAIST_URL,
            ),
            num_item(
                "whtr",
                "body",
                "腰高比 WHtR",
                waist / h,
                3,
                "",
                whtr_status,
                whtr_dir,
                "成人简便界值0.5",
                90,
                "https://pubmed.ncbi.nlm.nih.gov/20819243/",
            ),
        ]
    if waist is not None and hip is not None:
        items.append(
            num_item(
                "whr",
                "body",
                "腰臀比 WHR",
                whr_raw,
                3,
                "",
                whr_status,
                whr_dir,
                "男性低于0.90" if male else "女性低于0.85",
                100,
                "https://www.who.int/publications/i/item/9789241501491",
            )
        )
    if waist is not None:
        items += [
            num_item(
                "rfm",
                "body",
                "RFM估算体脂率",
                64 - 20 * h / waist + (0 if male else 12),
                1,
                "%",
                None,
                None,
                "身高和腰围估算",
                110,
                "https://pubmed.ncbi.nlm.nih.gov/30030479/",
            ),
        ]
    if hip is not None:
        items.append(
            num_item(
                "bai",
                "body",
                "BAI身体肥胖指数",
                hip / h_m**1.5 - 18,
                1,
                "",
                None,
                None,
                "个体误差可能较大",
                120,
                "https://pubmed.ncbi.nlm.nih.gov/26262696/",
            )
        )
    if waist is not None:
        items += [
            num_item(
                "absi",
                "body",
                "ABSI体型指数",
                (waist / 100) / (bmi_raw ** (2 / 3) * math.sqrt(h_m)),
                4,
                "",
                None,
                None,
                "原始值需结合人群标准化",
                130,
                "https://pubmed.ncbi.nlm.nih.gov/22815707/",
            ),
            num_item(
                "bri",
                "body",
                "BRI身体圆度指数",
                364.2
                - 365.5 * math.sqrt(1 - ((waist / 100 / (2 * math.pi)) ** 2) / (0.5 * h_m) ** 2),
                2,
                "",
                None,
                None,
                "无统一通用健康界值，用于体型比较",
                140,
                "https://pmc.ncbi.nlm.nih.gov/articles/PMC3692604/",
            ),
            num_item(
                "conicity",
                "body",
                "锥度指数",
                (waist / 100) / (0.109 * math.sqrt(w / h_m)),
                3,
                "",
                None,
                None,
                "无统一通用健康界值",
                150,
                "https://www.sciencedirect.com/science/article/pii/S2174204917301253",
            ),
        ]
    if not watson_flag:
        items += [
            num_item(
                "watson_tbw",
                "body",
                "Watson估算体水分",
                tbw,
                1,
                "L",
                None,
                None,
                "人体测量公式估算",
                160,
                WATSON_URL,
            ),
            num_item(
                "body_water_pct",
                "body",
                "估算体水分占比",
                tbw / w * 100,
                1,
                "%",
                None,
                None,
                "由Watson体水分除以体重推算",
                170,
                WATSON_URL,
            ),
        ]
    items.append(
        num_item(
            "blood_volume_nadler",
            "body",
            "Nadler估算血容量",
            blood,
            2,
            "L",
            None,
            None,
            "人体测量公式估算，不用于医疗决策",
            180,
            "https://pubmed.ncbi.nlm.nih.gov/21936146/",
        )
    )
    items += [
        num_item(
            "mifflin_ree",
            "metabolism",
            "Mifflin静息代谢",
            mifflin,
            0,
            "kcal/天",
            None,
            None,
            "静息状态估算",
            210,
            "https://pubmed.ncbi.nlm.nih.gov/2305711/",
        ),
        num_item(
            "harris_ree",
            "metabolism",
            "Harris–Benedict",
            harris,
            0,
            "kcal/天",
            None,
            None,
            "对照公式",
            220,
            "https://pubmed.ncbi.nlm.nih.gov/11817239/",
        ),
    ]
    if cunningham is not None:
        items.append(
            num_item(
                "cunningham_ree",
                "metabolism",
                "Cunningham",
                500 + 22 * (w - bf * w / 100),
                0,
                "kcal/天",
                None,
                None,
                "基于估算去脂体重",
                230,
                "https://pubmed.ncbi.nlm.nih.gov/7435418/",
                version="2.0",
            )
        )
    items.append(
        num_item(
            "tdee",
            "metabolism",
            "每日总能量 TDEE",
            mifflin * f["pal"],
            0,
            "kcal/天",
            f"PAL {f['pal']}",
            None,
            "没有统一健康范围",
            240,
            "https://www.fao.org/4/y5686e/y5686e07.htm",
        )
    )
    if rhr is not None:
        items.append(
            item(
                "resting_hr",
                "heart",
                "静息心率",
                str(rhr),
                rhr,
                "次/分钟",
                rhr_status,
                rhr_dir,
                "多数成人约60～100",
                310,
                AHA_HR_URL,
            )
        )
    items += [
        num_item(
            "max_hr_tanaka",
            "heart",
            "Tanaka最大心率",
            max_hr_raw,
            0,
            "次/分钟",
            None,
            None,
            "年龄公式估算",
            320,
            "https://pubmed.ncbi.nlm.nih.gov/11153730/",
        ),
        range_item(
            "moderate_hr",
            "heart",
            "中等强度心率",
            mod_lo,
            mod_hi,
            0,
            "次/分钟",
            "最大心率50%～70%",
            330,
            AHA_HR_URL,
        ),
        range_item(
            "vigorous_hr",
            "heart",
            "较高强度心率",
            vig_lo,
            vig_hi,
            0,
            "次/分钟",
            "最大心率70%～85%",
            340,
            AHA_HR_URL,
        ),
    ]
    if rhr is not None:
        items += [
            range_item(
                "karvonen_hr",
                "heart",
                "Karvonen目标心率",
                karv_lo,
                karv_hi,
                0,
                "次/分钟",
                "心率储备的50%～70%",
                350,
                CDC_KARVONEN_URL,
            ),
            num_item(
                "heart_rate_reserve",
                "heart",
                "心率储备",
                hrr_raw,
                0,
                "次/分钟",
                None,
                None,
                "估算最大心率减静息心率",
                360,
                CDC_KARVONEN_URL,
            ),
            num_item(
                "estimated_vo2max",
                "heart",
                "估算最大摄氧量",
                15.3 * max_hr_raw / rhr,
                1,
                "ml/kg/min",
                None,
                None,
                "使用估算最大心率，误差会叠加",
                370,
                "https://doi.org/10.1007/s00421-003-0988-y",
            ),
        ]
    items.append(
        range_item(
            "protein_range",
            "nutrition",
            "运动人群蛋白质参考",
            protein_lo,
            protein_hi,
            0,
            "g/天",
            "健康运动人群1.4～2.0 g/kg/天",
            410,
            "https://pubmed.ncbi.nlm.nih.gov/28642676/",
        )
    )
    if f["exercise_type"] is not None:
        items.append(
            item(
                "exercise_kcal",
                "exercise",
                "本次运动估算热量",
                str(kcal),
                kcal,
                "kcal",
                None,
                None,
                f"{label} × {f['exercise_minutes']}分钟",
                510,
                "https://pacompendium.com/adult-compendium/",
                metadata={
                    "exercise_type": f["exercise_type"],
                    "met": met,
                    "minutes": f["exercise_minutes"],
                },
                version="2.0",
            )
        )
    if f["bp_context"] is not None:
        items += [
            item(
                "systolic_bp",
                "vitals",
                "收缩压",
                str(sys_bp),
                sys_bp,
                "mmHg",
                sys_status,
                sys_dir,
                f"{ctx_label}；理想值低于120",
                610,
                NHC_BP_URL,
                metadata={"measurement_context": f["bp_context"]},
            ),
            item(
                "diastolic_bp",
                "vitals",
                "舒张压",
                str(dia_bp),
                dia_bp,
                "mmHg",
                dia_status,
                dia_dir,
                f"{ctx_label}；理想值低于80",
                620,
                NHC_BP_URL,
                metadata={"measurement_context": f["bp_context"]},
            ),
            item(
                "blood_pressure_category",
                "vitals",
                "血压判断",
                bp_cat,
                None,
                "",
                bp_cat,
                bp_dir,
                "单次读数仅供筛查，不能诊断",
                630,
                NHC_BP_URL,
                metadata={"measurement_context": f["bp_context"]},
            ),
        ]
    if f["sleep_hours"] is not None:
        items.append(
            item(
                "sleep_duration",
                "vitals",
                "平均每晚睡眠时长",
                str(f["sleep_hours"]),
                f["sleep_hours"],
                "小时",
                sleep_status,
                sleep_dir,
                sleep_ref,
                640,
                NHC_SLEEP_URL,
                metadata={"reference_low": sleep_lo, "reference_high": sleep_hi},
            )
        )

    tips, tips_en = build_tips(
        f,
        bmi_raw,
        bmi_status,
        wadj,
        waist_status,
        bf_flag,
        bf_out_of_age,
        watson_flag,
        rhr_status,
        bp_dir,
        sleep_dir,
    )

    return {
        "profile": {
            "age": age,
            "height_cm": h,
            "nickname": f["nickname"],
            "sex": sex,
            "weight_kg": w,
        },
        "body": {
            "absi": absi,
            "bai_pct": bai,
            "blood_volume_l": round(blood, 2),
            "bmi": bmi,
            "bmi_status": bmi_status,
            "body_fat_deurenberg_pct": round(deur, 1),
            "body_fat_deurenberg_status": deur_status,
            "body_fat_pct": round(bf, 1),
            "body_fat_status": bf_status,
            "body_water_pct": water_pct,
            "bri": bri,
            "bsa_m2": round(bsa, 2),
            "conicity": conicity,
            "fat_mass_kg": fat_mass,
            "ffm_kg": ffm,
            "janma_lbm_kg": round(janma, 1),
            "rfm_pct": rfm,
            "target_weight_high": tw_high,
            "target_weight_low": tw_low,
            "waist_status": waist_status,
            "watson_tbw_l": watson_tbw,
            "weight_adjustment_direction": (
                "gain" if wadj and wadj[1] else "lose" if wadj else None
            ),
            "weight_adjustment_kg": wadj[0] if wadj else None,
            "whr": whr,
            "whr_limit": whr_limit,
            "whr_status": whr_status,
            "whtr": whtr,
            "whtr_status": whtr_status,
        },
        "metabolism": {
            "cunningham_ree": cunningham,
            "harris_ree": round(harris, 0),
            "mifflin_ree": round(mifflin, 0),
            "pal": f["pal"],
            "tdee": tdee,
        },
        "heart": {
            "estimated_vo2max": vo2,
            "heart_rate_reserve": hrr,
            "karvonen_high": karv_hi,
            "karvonen_low": karv_lo,
            "max_hr": max_hr,
            "moderate_hr_high": mod_hi,
            "moderate_hr_low": mod_lo,
            "resting_hr": rhr,
            "resting_hr_status": rhr_status,
            "vigorous_hr_high": vig_hi,
            "vigorous_hr_low": vig_lo,
        },
        "nutrition": {"protein_high_g": protein_hi, "protein_low_g": protein_lo},
        "exercise": {
            "calories_kcal": kcal,
            "label": label,
            "met": met,
            "minutes": f["exercise_minutes"],
            "type": f["exercise_type"],
        },
        "vitals": {
            "blood_pressure_category": bp_cat,
            "blood_pressure_direction": bp_dir,
            "bp_context": f["bp_context"],
            "bp_context_label": ctx_label,
            "diastolic_bp": dia_bp,
            "diastolic_direction": dia_dir,
            "diastolic_status": dia_status,
            "sleep_direction": sleep_dir,
            "sleep_high": sleep_hi,
            "sleep_hours": f["sleep_hours"],
            "sleep_low": sleep_lo,
            "sleep_reference": sleep_ref,
            "sleep_status": sleep_status,
            "systolic_bp": sys_bp,
            "systolic_direction": sys_dir,
            "systolic_status": sys_status,
        },
        "items": items,
        "tips": tips,
        "tips_en": tips_en,
    }


# ---- record_id 持久计数器 -------------------------------------------------

DB_LOCK = threading.Lock()


def db_path():
    """计数器文件；HEALTH_DB 环境变量可覆盖（测试/多实例用）。"""
    return Path(os.environ.get("HEALTH_DB", "health-records.json"))


def next_record_id():
    """返回并递增持久化的 record_id（仅在提交成功时调用，失败不占号）。"""
    with DB_LOCK:
        path = db_path()
        try:
            n = json.loads(path.read_text())["next"]
        except (OSError, ValueError, KeyError, TypeError):
            n = 1
        with contextlib.suppress(OSError):
            path.write_text(json.dumps({"next": n + 1}) + "\n")
        return n
