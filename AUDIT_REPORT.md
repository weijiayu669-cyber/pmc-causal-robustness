# PMC 计算修正与数值审计报告（2026-11-21）

## Part A. Calculation verification

| Example | Old input | Corrected/verified input | Independent check |
|---|---|---|---|
| Education | partial_r=.29; rxz2=.224; ryz2=.177；错误交叉项为两个多重相关系数之积 | 矩阵交叉项 .18838965505016483；rxz2=.22352057677494624；ryz2=.17686271729239786；partial_r=.2984613894645406；nu=10562 | 完整4×4矩阵正定；逆精度矩阵及FWL残差偏相关一致；T(0,0)=32.138154101648，四舍五入32.1 |
| Retention 正式工作输入 | partial_r=−.13；rxz2=.2704；ryz2=.2025；nu=7635 | 保留；共享inputs.R，明确为假设依赖的四舍五入工作输入 | T(0,0)=−11.456422589690；与−13.27绝对差1.813577410310；不能声称恢复原文多层模型t |

Education 用 `rXZ' solve(RZZ) rYZ` 代替错误交叉项；使用完整精度参与计算，正文可展示0.30。完整矩阵的正定性、精度矩阵公式、实际回归残差和直接OLS均有测试。

## Part B. Updated numerical results

所有主结果均为100次独立完整Monte Carlo repetition的条件比例均值，每个参数1,000 draws。MCSE=跨repetition SD/10。下表不混入其他Retention重建情景。

| Example | Prior | Old probability | New probability | MCSE | Mean admissible |
|---|---|---:|---:|---:|---:|
| Retention | U(-1,1) | 0.678885694346454 | 0.678885694346 | 0.000872590642718 | 594284.2 |
| Retention | U(-0.5,0.5) | 0.83659422 | 0.83659422 | 0.000652365187216 | 1000000 |
| Education | U(-1,1) | 0.820084 | 0.825283992413 | 0.000781046584493 | 599629.39 |
| Education | U(-0.5,0.5) | 0.995132 | 0.99697109 | 3.97324285412e-05 | 1000000 |

跨100次运行SD依次为：0.00872590642718, 0.00652365187216, 0.00781046584493, 0.000397324285412。Education相对于用户给出的旧展示值变化为约+0.005199992413386与+0.00183909；旧CSV保留在审计工作记录中。Figure 2的标签为0.825、0.997，后者小于1。

A.6正式工作输入没有改变，仍完整重跑两次并核验一致：

| Prior | NHST Q0 | Effect-size Q_delta | Q_delta MCSE | Mean admissible |
|---|---:|---:|---:|---:|
| U(-1,1) | 0.678885694346 | 0.647879036366 | 0.000885235356829 | 594284.2 |
| U(-0.5,0.5) | 0.83659422 | 0.79485314 | 0.000694559479316 | 1000000 |

## Part C. Retention diagnostics

源表逐格核对：29个非空cell的Table 8人数、logit均值/SD，以及Table 9 reading均值/SD全部匹配原`mcdat()`。Table 9有9个cell人数不同，总差31人：retained L9 46→45、L14 39→38；promoted L1 3054→3044、L2 1661→1654、L3 983→977、L4 323→321、L5 441→440、L6 154→153、L7 213→211。Table 8合计471+7168=7639，Table 9合计469+7139=7608。空的promoted L15按缺失均值/SD保存，没有填0。

以下每行用1,000个合成数据集，seed123，每情景重置，共用同一计数设计内的潜在正态随机数。T按每次partial_r非线性计算后求均值，不是将平均r直接代入。nu采用PMC增广模型N−4；仅拟合y~x+z的普通回归df是N−3。这里的nu都不是原文多层模型df的恢复。

| Count reconstruction | rho_within | mean rxz | mean ryz | mean rxy | mean partial r | mean T(0,0) | abs(mean T+13.27) | mean abs(T+13.27) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| reference | -0.5 | 0.523155244 | -0.523488194 | -0.332686109 | -0.081009481 | -7.103343727 | 6.166656273 | 6.166656273 |
| reference | -0.25 | 0.523155244 | -0.487080951 | -0.332707565 | -0.104647909 | -9.196391834 | 4.073608166 | 4.073608166 |
| reference | 0 | 0.523155244 | -0.450594092 | -0.332697099 | -0.127450205 | -11.230310456 | 2.039689544 | 2.054398984 |
| reference | 0.25 | 0.523155244 | -0.414034221 | -0.332653452 | -0.149592518 | -13.222674437 | 0.047325563 | 0.845648586 |
| reference | 0.5 | 0.523155244 | -0.377401777 | -0.332573804 | -0.171223044 | -15.188650512 | 1.918650512 | 1.938741689 |
| outcome_matched | -0.5 | 0.522460241 | -0.523091882 | -0.332130959 | -0.080962817 | -7.084861620 | 6.185138380 | 6.185138380 |
| outcome_matched | -0.25 | 0.522460241 | -0.486593139 | -0.332255157 | -0.104750353 | -9.186848296 | 4.083151704 | 4.083151704 |
| outcome_matched | 0 | 0.522460241 | -0.450076619 | -0.332377166 | -0.127695741 | -11.229518400 | 2.040481600 | 2.066367485 |
| outcome_matched | 0.25 | 0.522460241 | -0.413545282 | -0.332495734 | -0.149979238 | -13.230825464 | 0.039174536 | 0.846921188 |
| outcome_matched | 0.5 | 0.522460241 | -0.377000741 | -0.332610572 | -0.171752595 | -15.206269532 | 1.936269532 | 1.961199608 |

正式四舍五入输入另列：r_xz=.52，r_yz=−.45，r_xy=−.333163331932726，partial_r=−.13，T0=−11.456422589690，绝对差1.813577410310。它与rho=0重建均值有舍入差异。

所有上述重建都**不能识别真实偏相关或原文t**。源表能识别各cell人数和边际均值/SD（受发表舍入限制），不能识别cell内Z–Y协方差、缺失阅读者的联合分布、阅读观测子样本logit矩、学校标识或测评间隔。参考重建把Table 9边际矩外推给缺失的31人；outcome-matched重建则把Table 8 logit矩沿用到阅读有观测的子样本。两者都需要假设。正态生成保持的是分布参数，并不强制每个有限样本恰等于发表矩。

rho网格在分析前固定，共同rho只是敏感性设定。rho=.25的均值接近−13.27不证明该rho真实，也不建立模型等价性。本次没有反推rho、partial_r或通过拟合原文t选择正式输入。计数改变的模拟差异约为0.00025 partial r，但在1,000次模拟MCSE约0.00037的量级，不能把这个有限模拟差异当成精确计数效应。CSV还提供不含模拟噪声的混合分布矩：rho=0时reference partial_r=−.127609178045541、T0=−11.2422072435671，matched为−.127837467417113、−11.2397652413535；这些也只是所设混合分布的结果。

原文Equation 4包含14个共享分层指标、连续propensity logit、测评间隔，以及相关的学校随机截距和retention随机斜率。仅控制一个z的普通回归近似与其不同。1–29的生成器strata只是cell编码，不能作为原文15个共享分层的控制变量使用。原文−9.01/.68因展示舍入不精确等于−13.27，不应用两位小数的商覆盖其报告t。原始论文见[DOI](https://doi.org/10.3102/01623737027003205)及[完整原文文本镜像](https://www.academia.edu/29274731/Effects_of_Kindergarten_Retention_Policy_on_Childrens_Cognitive_Growth_in_Reading_and_Mathematics)。源表核查基于文本镜像，未取得可视觉核对的原版PDF；详见SOURCE_AUDIT.md。

## Part D. Manuscript update checklist

仓库不含论文源文件，以下为需要应用到稿件的明确清单，不声称已经改过正文。保持现有符号与PMC框架。

- **正文Education**：把偏相关计算的交叉项改为r_XZᵀR_ZZ⁻¹r_YZ，两个平方多重相关分别用对应二次型；输入由0.29改为0.2984613894645406（展示0.30），0.224/0.177改用未舍入矩阵结果。写明T0=32.138154101648，与32.1的展示值相符。
- **正文Education结果及Abstract**：若出现旧0.820084、0.995132，分别替换为0.825283992413386、0.99697109；可写82.5284%和99.6971%，不得写严格100%。若Abstract仅用两位小数概率，教育restricted也应至少保留0.997。Retention正式概率保持0.678885694346454、0.83659422，删去其零混杂t等于−13.27的说法。
- **Appendix A.2**：若该处给出偏相关推导/计算，替换错误多控制交叉项，加入本报告Part A的四个矩阵验证量。保持T(0,0)=r_partial sqrt(nu/(1−r_partial²))及既有增广nu约定，不能把该nu称作源文mixed-model df。
- **Appendix A.4**：若该处描述Retention重建，逐格引用retention_source_cells.csv；区分Table 8 N7639与Table 9 N7608、31个缺失reading和两种外推假设；说明rho=0为参考假设，加入Part C固定网格诊断。说明一个连续z的OLS不等价于原文Equation 4。
- **Appendix A.5 / Table A3**：替换Education两行的概率、SD、MCSE、mean/min/max admissible与半程统计，完整精度来自robustness_summary.csv；前两行Retention经重跑不变。四行主概率、MCSE、mean admissible见Part B。方法说明MCSE=SD/sqrt(100)，Cartesian pairs不是独立Bernoulli样本。
- **Appendix A.6 / Table A4**：输入和结果不变，经两次完整重跑确认；NHST Q0、effect-size Q_delta和MCSE见Part B。保留Delta_min=.20×13.53=2.706、d_delta=Delta_min/SE_f、p_delta=Phi(−tt−d_delta)，不额外减1.96。SE_f=(13.53/.240548)sqrt(D/7635)/vx。注明13.53来自Table 2更大reading样本N11797，是工作尺度，不是已知at-risk或mixed-model残差SD。
- **Figure 1**：以重跑输出替换PDF/PNG；曲面正式输入未变，标签更新为0.679、0.837，并清楚标注假设依赖。
- **Figure 2**：以修正输入重绘曲面及V，PDF/PNG标签0.825、0.997；保留黑白风格及现有图例。

A.2/A.4具体段落标题未随仓库提供；以上按所含公式/重建内容定位，不捏造稿件原句、公式号或页码。其余PMC六步、符号、100 repetitions、1,000 draws、两组uniform scenarios与seed123均保留。

## Part E. Changed files and validation

修改：robustness.R、minimum_effect.R、tests/minimum_effect.R、README.md、robustness_summary.csv、education.pdf/png、retention.pdf/png。A.6两个CSV重跑后与原版数值和文件字节一致，故Git不产生内容diff。

新增：inputs.R；retention_reconstruction.R；retention_audit.R；robustness_runs.csv；retention_source_cells.csv；retention_cell_source_audit.csv；retention_diagnostics.csv；tests/calculation_audit.R；tests/robustness_results.R；tests/retention_audit.R；两份alignment-layout.json；SOURCE_AUDIT.md；本报告；reproducibility.txt。

实际环境R 4.5.0 (2025-04-11)，aarch64-apple-darwin20，macOS27.0.1。分析只用base/recommended包，版本与R一致；seed123，Mersenne-Twister/Inversion/Rejection。执行命令和完整sessionInfo见reproducibility.txt及README。

验证：

1. 修复前的Education .29对精确偏相关断言失败；修复后calculation_audit通过，包括4×4正定、精度矩阵、FWL与直接OLS。
2. retention_audit测试通过；源表全部边际及人数检查、缺失outcome、默认生成器逐位相同与RNG状态相同；高样本量敏感性边际/相关及混合矩独立验证。
3. 两次完整robustness运行CSV逐字节一致；独立重放全部400次repetition，所有合法pair的T/probability有限、V一致，超阈值计数与每次比例相符；100-run均值、SD和MCSE检查通过。
4. minimum_effect点对点OLS、合法区域和阈值测试通过；两次200-row运行CSV逐字节一致，独立FWL和效果尺度事件检查全部200个NHST/effect-threshold结果通过。
5. 两次完整2×5×1000个Retention合成数据集审计CSV逐字节一致；这是输入诊断，不混进100-repetition主分析。
6. PDF两图面板对齐通过1.5pt容差，最小字9pt，碰撞审计0 fail/0 warn；PNG目检保留原风格。通用Nature静态检查的SVG/TIFF/Cairo偏好与本任务保留PDF/PNG的要求不同，未据此添加格式或改设计；不宣称符合某期刊全部投稿规范。

尚未解决的是源表不识别联合分布、缺失机制与mixed-model所需数据，及原PDF视觉核验的访问限制。本次没有恢复原研究完整个体数据，也未改变正式Retention输入以追求−13.27。
