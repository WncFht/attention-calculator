# 第三方组件说明

`scripts/fetch-vendor-assets.sh` 拉取的前端依赖（存于已 gitignore 的 `src/attention_calculator/static/vendor/`）：

| 组件                                      | 来源                    | 许可证     |
| ----------------------------------------- | ----------------------- | ---------- |
| MathJax 3（tex-mml-chtml.js、tex-svg.js） | jsdelivr `mathjax@3`    | Apache-2.0 |
| KaTeX 0.16.9（css、js、woff2 字体）       | jsdelivr `katex@0.16.9` | MIT        |
| html2canvas                               | html2canvas.hertzen.com | MIT        |

`src/attention_calculator/static/` 下的 `bg1.png`、`bg2.png`、`title.png` 与 `templates/` 下各页面为 zhuyidao.net 站端资产的逐字节镜像，仅用于复现研究，版权归原作者所有。
