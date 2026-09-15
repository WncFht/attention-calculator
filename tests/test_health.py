"""Tests for the /health sibling app clone (GET /health*, POST /health/calculate).

行为全部以 bench/data/health-probes.jsonl 的线上实测为准（byte-exact 由
bench/parity_health.py 负责）；这里钉死契约骨架：路由面、包络形状、校验
顺序与错误文案、公式/阈值关键点、record_id 只认成功提交。
"""

import json

import pytest

from attention_calculator import server

REQ = {
    "nickname": "Probe",
    "sex": "male",
    "age": 35,
    "height_cm": 175,
    "weight_kg": 70,
    "pal": 1.55,
}


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("HEALTH_DB", str(tmp_path / "records.json"))
    return server.app.test_client()


def calc(client, **over):
    """POST /health/calculate with REQ updated by ``over``; None removes."""
    payload = dict(REQ)
    for k, v in over.items():
        if v is None:
            payload.pop(k, None)
        else:
            payload[k] = v
    return client.post("/health/calculate", json=payload)


def test_pages_and_route_edges(client):
    """路由面实测：/health 与 /health/ 同页面；/health/en 独立英文模板；
    /health/en/ 与域内未匹配路径一律 500；/healthxyz 落主站 404。"""
    zh = client.get("/health").get_data()
    assert zh == client.get("/health/").get_data()
    assert len(zh) > 10000 and "健康" in zh.decode()
    en = client.get("/health/en").get_data()
    assert len(en) > 10000 and en != zh
    for path in ("/health/en/", "/health/xyz", "/health/static/bg1.png"):
        resp = client.get(path)
        assert resp.status_code == 500, path
        assert resp.get_json() == {"error": "服务器内部错误，请稍后再试。", "ok": False}
    resp = client.get("/healthxyz")
    assert resp.status_code == 404
    assert resp.get_json() == {"error": "请求的页面不存在"}  # 主站裸包络
    assert client.get("/health/calculate").status_code == 500
    assert client.open("/health/calculate", method="OPTIONS").status_code == 200
    assert client.post("/health/").status_code == 500
    assert client.post("/health/en").status_code == 500


def test_transport_errors(client):
    """非 JSON / 坏 JSON / 非 dict 顶层 -> 400 同一文案（tr:* 探针）。"""
    for kw in (
        {"data": "nickname=x"},  # form CT
        {"data": "{", "content_type": "application/json"},
        {"data": "null", "content_type": "application/json"},
        {"data": "[]", "content_type": "application/json"},
        {"data": '"hi"', "content_type": "application/json"},
    ):
        resp = client.post("/health/calculate", **kw)
        assert resp.status_code == 400, kw
        assert resp.get_json() == {"error": "提交内容格式不正确。", "ok": False}
    # charset 后缀的 JSON CT 仍算 JSON（s2:ct-charset）
    resp = client.post(
        "/health/calculate", data=json.dumps(REQ), content_type="application/json; charset=utf-8"
    )
    assert resp.status_code == 200


def test_validation_order_and_messages(client):
    """错误文案逐字 + 校验顺序（成对探针钉死）。"""
    cases = [
        ({"nickname": None}, "请填写昵称。"),
        ({"nickname": "a" * 25}, "昵称不能超过24个字符。"),  # 长度先于字符集
        ({"nickname": "a" * 24 + "@"}, "昵称不能超过24个字符。"),
        ({"nickname": "a@b"}, "昵称只能包含中文、字母、数字、下划线或短横线。"),
        ({"sex": None}, "请选择性别，以便使用对应公式。"),
        ({"sex": "other"}, "请选择性别，以便使用对应公式。"),
        ({"sex": True}, "请选择性别，以便使用对应公式。"),
        ({"age": None}, "请填写年龄。"),
        ({"age": "abc"}, "年龄必须是数字。"),
        ({"age": True}, "年龄格式不正确。"),
        ({"age": 17}, "年龄应在 18～120 之间。"),
        ({"age": 17.9}, "年龄应在 18～120 之间。"),  # 范围先于整数
        ({"age": 119.5}, "年龄必须是整数。"),
        ({"height_cm": None}, "请填写身高。"),
        ({"height_cm": "abc"}, "身高必须是数字。"),
        ({"height_cm": True}, "身高格式不正确。"),
        ({"height_cm": 0.9}, "身高应在 1～250 之间。"),
        ({"weight_kg": 0.9}, "体重应在 1～500 之间。"),
        ({"pal": None}, "请填写活动水平。"),
        ({"pal": True}, "活动水平格式不正确。"),
        ({"pal": 1.3}, "活动水平应在 1.4～2.4 之间。"),
        ({"pal": 1.4}, "请选择页面提供的活动水平。"),  # 范围内但非三档枚举
        ({"waist_cm": 0.9}, "腰围应在 1～300 之间。"),
        ({"waist_cm": True}, "腰围格式不正确。"),
        ({"waist_cm": [80]}, "腰围必须是数字。"),
        ({"hip_cm": 300.1}, "臀围应在 1～300 之间。"),
        ({"resting_hr": "abc"}, "静息心率必须是数字。"),
        ({"resting_hr": 19.5}, "静息心率应在 20～250 之间。"),  # 范围先于整数
        ({"resting_hr": 62.5}, "静息心率必须是整数。"),
        ({"exercise_minutes": 30}, "运动类型和运动时长需要一起填写，或都留空。"),
        ({"exercise_type": "bogus"}, "请选择页面提供的运动类型。"),  # 枚举先于成对
        ({"exercise_type": "jogging", "exercise_minutes": 0.5}, "运动时长应在 1～1440 之间。"),
        ({"exercise_type": "jogging", "exercise_minutes": 30.5}, "运动时长必须是整数。"),
        ({"exercise_type": "jogging", "exercise_minutes": True}, "运动时长格式不正确。"),
        ({"bp_context": "bogus"}, "请选择页面提供的血压测量场景。"),  # 枚举先于成组
        ({"systolic_bp": 118}, "血压测量场景、收缩压和舒张压需要一起填写，或都留空。"),
        (
            {"bp_context": "clinic", "systolic_bp": 49.5, "diastolic_bp": 76},
            "收缩压应在 50～300 之间。",
        ),
        ({"bp_context": "clinic", "systolic_bp": 118.7, "diastolic_bp": 76}, "收缩压必须是整数。"),
        (
            {"bp_context": "clinic", "systolic_bp": 85, "diastolic_bp": 95},
            "收缩压通常应高于舒张压，请检查输入。",
        ),
        (
            {"bp_context": "clinic", "systolic_bp": 100, "diastolic_bp": 100},
            "收缩压通常应高于舒张压，请检查输入。",
        ),  # <= 拒
        ({"sleep_hours": "abc"}, "平均每晚睡眠时长必须是数字。"),
        ({"sleep_hours": 24.1}, "平均每晚睡眠时长应在 0～24 之间。"),
        # 跨字段顺序：腰围先于运动成对，运动先于血压成组
        ({"waist_cm": 0, "exercise_type": "jogging"}, "腰围应在 1～300 之间。"),
        (
            {"exercise_type": "jogging", "bp_context": "clinic"},
            "运动类型和运动时长需要一起填写，或都留空。",
        ),
    ]
    for over, err in cases:
        resp = calc(client, **over)
        assert resp.status_code == 400, over
        assert resp.get_json() == {"error": err, "ok": False}, over


def test_honeypot(client):
    """website 蜜罐：真值即拒，且先于一切字段错误（ord:* 探针）。"""
    resp = calc(client, website="x", nickname=None)
    assert resp.get_json() == {"error": "提交未通过校验。", "ok": False}
    assert calc(client, website="").status_code == 200  # 空串放行
    assert calc(client, website=None).status_code == 200  # null 同缺席


def test_success_shape_and_record_id(client):
    """成功包络 {ok,record_id,results}；record_id 持久递增、失败不占号。"""
    resp = calc(client)
    assert resp.status_code == 200
    body = resp.get_json()
    assert body["ok"] is True and isinstance(body["record_id"], int)
    res = body["results"]
    assert set(res) == {
        "profile",
        "body",
        "metabolism",
        "heart",
        "nutrition",
        "exercise",
        "vitals",
        "items",
        "tips",
        "tips_en",
    }
    assert len(res["items"]) == 19  # 仅必填 -> 19 张卡片
    assert res["profile"] == {
        "age": 35,
        "height_cm": 175.0,
        "nickname": "Probe",
        "sex": "male",
        "weight_kg": 70.0,
    }
    first = body["record_id"]
    assert calc(client).get_json()["record_id"] == first + 1  # 成功 +1
    calc(client, age=17)  # 失败不占号
    assert calc(client).get_json()["record_id"] == first + 2


def test_core_formulas(client):
    """公式/状态抽查：bmi(raw 判状态)、wadj、mifflin、tdee(raw 底)、whtr。"""
    res = calc(client, height_cm=200, weight_kg=73.9).get_json()["results"]
    assert res["body"]["bmi"] == 18.5  # 显示舍入
    assert res["body"]["bmi_status"] == "体重过低"  # 状态看 raw 18.475
    assert res["body"]["weight_adjustment_kg"] == 0.1  # 18.5*4-73.9
    assert res["body"]["weight_adjustment_direction"] == "gain"

    res = calc(client, height_cm=170, weight_kg=70, pal=1.85).get_json()["results"]
    m = res["metabolism"]
    assert m["mifflin_ree"] == 1592.0  # raw 1592.5 banker's
    assert m["tdee"] == 2946.0  # raw*1.85=2946.125
    assert m["harris_ree"] == 1643.0  # raw 1643.287
    assert m["cunningham_ree"] == 1708.0  # 500+22*ffm_raw 54.9

    res = calc(client, waist_cm=87.5).get_json()["results"]  # 87.5/175 = 0.5
    assert res["body"]["whtr"] == 0.5
    assert res["body"]["whtr_status"] == "达到或超过0.5参考界值"


def test_bf_flag_and_cascade(client):
    """CUN-BAE/Deurenberg 不宜界值 [0,70]；CUN-BAE 旗标停算派生项。"""
    res = calc(client, height_cm=200, weight_kg=27).get_json()["results"]
    assert res["body"]["body_fat_status"] == "结果超出公式的合理解释范围"
    assert res["body"]["fat_mass_kg"] is None  # 级联停算
    assert res["metabolism"]["cunningham_ree"] is None
    assert "CUN-BAE给出了不宜解释的结果" in res["tips"][2]

    # 超龄只换状态文案，不停派生计算（a81 实测 fm/cunn 仍出）
    res = calc(client, age=81).get_json()["results"]
    assert res["body"]["body_fat_status"] == "超出CUN-BAE原始18～80岁验证范围"
    assert res["body"]["fat_mass_kg"] is not None
    assert "CUN-BAE原始研究对象为18～80岁" in res["tips"][-1]

    # Deurenberg 单侧旗标不影响 CUN-BAE 派生
    res = calc(client, height_cm=170, weight_kg=200).get_json()["results"]
    assert res["body"]["body_fat_deurenberg_status"] == "结果超出公式的合理解释范围"
    assert res["body"]["body_fat_status"] == "达到高体脂报警界值"


def test_watson_flag(client):
    """Watson 体水分 >体重（pct>100）-> 两字段缺席 + 提示。"""
    res = calc(client, age=80, height_cm=250, weight_kg=30, sex="female").get_json()["results"]
    assert res["body"]["watson_tbw_l"] is None
    assert res["body"]["body_water_pct"] is None
    assert any("Watson" in t for t in res["tips"])
    keys = {i["key"] for i in res["items"]}
    assert "watson_tbw" not in keys and "body_water_pct" not in keys


def test_bp_categories(client):
    """血压类别高侧优先；家庭场景无分级。"""

    def cat(s, d, ctx="clinic"):
        return calc(client, bp_context=ctx, systolic_bp=s, diastolic_bp=d).get_json()["results"][
            "vitals"
        ]

    assert cat(118, 76)["blood_pressure_category"] == "正常血压范围"
    assert cat(120, 80)["blood_pressure_category"] == "正常高值范围"
    assert cat(150, 50)["blood_pressure_category"] == "1级高血压范围"  # 高侧压偏低
    assert cat(160, 70)["blood_pressure_category"] == "2级高血压范围"
    assert cat(180, 50)["blood_pressure_category"] == "3级高血压范围"
    assert cat(95, 35)["blood_pressure_category"] == "血压偏低范围"
    assert cat(160, 100, "home")["blood_pressure_category"] == "达到家庭血压高血压筛查界值"
    assert cat(140, 50, "home")["blood_pressure_category"] == "达到家庭血压高血压筛查界值"
    assert cat(118, 76, "home")["bp_context_label"] == "家庭自测"


def test_sleep_and_hr(client):
    """睡眠参考按年龄切换；心率区间与 vo2max 抽查。"""
    res = calc(client, sleep_hours=6.5).get_json()["results"]
    assert res["vitals"]["sleep_status"] == "低于睡眠时长参考范围"
    assert res["vitals"]["sleep_reference"] == "成年人一般7～8小时"
    res = calc(client, age=65, sleep_hours=6.5).get_json()["results"]
    assert res["vitals"]["sleep_status"] == "睡眠时长参考范围"
    assert res["vitals"]["sleep_reference"] == "老年人一般6～7小时"

    res = calc(client, resting_hr=250).get_json()["results"]
    h = res["heart"]
    assert h["max_hr"] == 184.0  # round(183.5) banker's
    assert h["heart_rate_reserve"] == -66.0
    # Karvonen 用 raw 最大心率且不交换（站端怪癖：低>高）
    it = {i["key"]: i for i in res["items"]}
    assert it["karvonen_hr"]["display_value"] == "217.0～203.0"


def test_items_shape(client):
    """卡片键集与 sort_order 契约（sibling-apps.md §2）。"""
    res = calc(
        client,
        waist_cm=82,
        hip_cm=96,
        resting_hr=62,
        exercise_type="jogging",
        exercise_minutes=30,
        bp_context="clinic",
        systolic_bp=118,
        diastolic_bp=76,
        sleep_hours=7.5,
    ).get_json()["results"]
    orders = [i["sort_order"] for i in res["items"]]
    assert orders == sorted(orders)  # 卡片按 sort_order 序
    assert {i["key"] for i in res["items"]} == {
        "bmi",
        "target_weight",
        "body_fat_cun_bae",
        "body_fat_deurenberg",
        "fat_mass",
        "ffm",
        "janma_lbm",
        "bsa_mosteller",
        "waist_assessment",
        "whtr",
        "whr",
        "rfm",
        "bai",
        "absi",
        "bri",
        "conicity",
        "watson_tbw",
        "body_water_pct",
        "blood_volume_nadler",
        "mifflin_ree",
        "harris_ree",
        "cunningham_ree",
        "tdee",
        "resting_hr",
        "max_hr_tanaka",
        "moderate_hr",
        "vigorous_hr",
        "karvonen_hr",
        "heart_rate_reserve",
        "estimated_vo2max",
        "protein_range",
        "exercise_kcal",
        "systolic_bp",
        "diastolic_bp",
        "blood_pressure_category",
        "sleep_duration",
    }
    for i in res["items"]:
        assert set(i) == {
            "key",
            "category_key",
            "sort_order",
            "display_name",
            "display_value",
            "numeric_value",
            "unit",
            "status",
            "direction",
            "reference",
            "source_url",
            "formula_version",
            "metadata",
        }
