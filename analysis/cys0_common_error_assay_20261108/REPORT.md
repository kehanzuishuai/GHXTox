# 为什么 Cys=0 强漏判会跨模型重复

以六种子纯 ProtT5 的已冻结定义取 67 条原标签为有毒、Cys=0、校准后有毒概率不超过 0.10 的训练集留出样本。六个种子的判断方向均为无毒。复用现有 EGNN、纯 ESM-2、历史 PepBERT 两个规格、ESM-IF1 线性探针、ESM-2 加 ESM-IF1 六个分支，在 0.5 参考点上有 62 条全部判为无毒。各分支训练方案和校准方式不同，这一数字用于定义本次证据核查范围。未使用 Test1/Test2，也未训练或改标签。

## 逐条核查结果

| 首要证据类型 | 67 条重点样本 | 62 条共同漏判 |
|---|---:|---:|
| 毒蛋白或完整毒素的部分序列 | 22 | 21 |
| 输入为前体，原文鉴定或测量的是成熟肽 | 4 | 4 |
| 数据库读数指向不同序列的类似物 | 1 | 1 |
| 对应分子有直接宿主红细胞或细胞损伤证据 | 10 | 8 |
| 特定功能或条件效应 | 15 | 13 |
| 身份有线索，直接实验未核实 | 6 | 6 |
| 仅找到条件受限的弱或阴性结果 | 1 | 1 |
| 未确认分子身份 | 8 | 8 |

类别互斥，只用于组织证据。有直接宿主损伤的样本仍可能存在修饰不明或测量单位待核实。实验档案按精确裸序列匹配到 13 条目标、52 条记录，其中数据库回链、同序列的不同化学形式及重复测量均保留原级别。CAPTP 的公开训练文件仅有序列与二元标签，逐条缺少从该标签回到具体实验的明确记录。外部原文可支持相应实验结果，不能替代原标签谱系。 [CAPTP 数据说明](https://academic.oup.com/bioinformatics/article/40/5/btae297/7663469) · [旧文件一致性校验](../ghxtox_error_atlas_20261025/evidence_phase2/captp_source_verification.json)

### 上游标签能追到哪里

逐条对齐 CAPTP 声明合并的三个上游来源。67 条中有 60 条精确匹配至少一份上游公开清单。CSM-Toxin 的 UniProt 阳性清单覆盖 54 条主清单、4 条盲测清单，盲测清单的这 4 条也在主清单中。进一步按其 CD-HIT 70% 文件核对，主清单命中 44 条，盲测清单命中 4 条。ATSE 原训练集和原测试集共命中 33 条，全部原标阳性。ToxinPred2 主阳性集命中 4 条。各来源会重叠。尤其是 31 aa 的 BmLAO 记录和 21 aa 的 Granulitoxin 片段，在 CSM-Toxin 或 ATSE 的阳性清单里已经存在。这个结果能说明对象不一致早于当前模型出现。无法由名单重合确定 CAPTP 合并时采用了哪一份记录，或当时对应哪个实验。 [CSM-Toxin 数据与代码](https://bitbucket.org/ascherslab/csm-toxin/src/c5daaa421a2bb21d902dd3b64a2f24bb2e7f0877/data_processing/) · [ATSE 数据](https://github.com/WLYLab/ATSE/tree/main/Data) · [ToxinPred2 数据](https://webs.iiitd.edu.in/raghava/toxinpred2/stand.html)

### 分子对象不一致有直接证据

- **Granulitoxin，peptide_162。** 输入是完整分离物的 N 端 21 aa。小鼠腹腔给药的 2 h LD50 400 ± 83 μg/kg 测的是 4,958 Da 完整分离物。类似的例子包括 pleuturegin 的 28 aa 片段和 Aralin A 链的 25 aa 片段。 [Granulitoxin 原文](https://pubmed.ncbi.nlm.nih.gov/9876306/) · [Pleuturegin 原文](https://pubmed.ncbi.nlm.nih.gov/11676502/) · [Aralin 原文](https://pubmed.ncbi.nlm.nih.gov/12963119/)
- **BmLAO，peptide_11。** 输入的 31 aa 恰好拼接了原文分别合成并各自 C 端酰胺化的 13、8、10 aa 三条衍生物。论文列出的红细胞实验对象是各条衍生物或完整蛋白，旧数据库对 31 aa 的读数回链不成立。该样本六种子判无毒，但 PepBERT-large 在 0.5 判有毒，计在 67 条而非 62 条共同漏判。 [原论文 Results 与 Methods](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0033639)
- **Ponericin-W3，peptide_259。** 原文详细实验使用末端少一个 Lys 的 W3-desK，旧数据库将马红细胞的结果关联到天然 W3。该实验支持类似物的结果，不能直接给天然 W3 加上同一实测值。 [Orivel 原论文](https://doi.org/10.1074/jbc.M100216200)
- **Dv63、Dv73、Dv5 与 Mri1a。** 四条输入带前体区域，论文测的是较短成熟体。Dv63 的果蝇麻痹 EP50 为 10.27 nmol/g，Dv73 的线虫 EC50 为 94.5 μM，Dv5 在多个所测终点为弱或未检出。Mri1a 蚜虫实验测的是酰胺化十肽。应分别保留分子边界和测量终点。 [Dv 肽原论文 Table 2](https://doi.org/10.1073/pnas.2023815118) · [Mri1a 原论文](https://pubmed.ncbi.nlm.nih.gov/31557881/)

### 共同漏判里仍有可信的直接实验

- **同一人红细胞实验。** 三条均被所有分支判为无毒。酰胺化 EMP-AF、Mastoparan-L、Polybia-MP III 在 320 μg/mL、PBS、37°C、1 h 下分别为 100%、100%、79% 溶血。该实验中，GHXTox 原无毒的 Protonectarina-MP、Mastoparan-J、PMM2 对应实验分子也均为 100%。相同测量条件下，这些结果说明原广义标签与实测溶血可不一致。 [实验原文 Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10611374/) · [旧面板与标签匹配](../cys0_membrane_evidence_20261107/REPORT.md)
- **化学形式在实验里可改变效应。** 同序列 Mastoparan-L 末端酰胺化比游离酸有更强细胞毒性。Conomap-Vt 的第 2 位 D-Phe 增强神经兴奋相关活性。Contulakin-G 同骨架的糖基化版本在小鼠脑室给药中约低 10 倍剂量就出现相应运动效应。这些是输入确实缺信息的证据，尚未测得补入后能纠正多少 GHXTox 错题。 [Mastoparan-L Fig. 4A](https://pmc.ncbi.nlm.nih.gov/articles/PMC5097029/) · [Conomap-Vt 原文](https://doi.org/10.1016/j.febslet.2006.06.011) · [Contulakin-G 原文](https://pubmed.ncbi.nlm.nih.gov/10318778/)
- **毒性的含义随终点而变。** Nc2a 的 HEK293 CC50 为 48.2–56.9 μM，人红细胞 HC50 >250 μM。Polybine-II 有炎症效应，原研究条件下未见溶血。CPY-Pl1 阻断 Kv1.6，IC50 170 nM。Mr-1 原文在 100 μM 测四种大鼠 nAChR 亚型，仅见有限读数，作者总结无显著活性。这些条件各自独立。 [Nc2a 原文](https://pubmed.ncbi.nlm.nih.gov/34302796/) · [Polybine-II 原文](https://pubmed.ncbi.nlm.nih.gov/15572194/) · [CPY-Pl1 Fig. 6](https://pubmed.ncbi.nlm.nih.gov/18505731/) · [Mr-1 原文](https://www.scielo.br/j/jvatitd/a/6Xtt7D6swjH7h4XGj33LwYS/?format=pdf)

### 数据分布能解释什么

在 6,387 条训练肽中，Cys=0 的原有毒比例为 207/3,213（6.4%），Cys>0 则为 1,611/3,174（50.8%）。原有毒样本中强反向比例为 67/207（32.4%）和 21/1,611（1.3%），五折都出现同方向差异。去掉全部旧片段标记后，仍为 45/141（31.9%）和 12/1,359（0.9%）。这说明已知片段不足以解释整体关联。无片段标记组仍包括本次新查出的前体和类似物，不能称为核实过的干净标签。家族、来源和分子形式也未独立控制。之前同折、同标签、同 Cys、相近长度的匹配只有 43/67 条有参照，pLDDT 和 Cα 接触密度差的方向跨折不稳定。现有结果支持把数据分布、实验终点和输入分子分别验证，无法把模型行为归因于单一因素。 [六种子分型与匹配](../six_seed_1d_weakness_20261107/REPORT.md)

### 下一步的实验顺序

1. 先用有明确实测分子的案例和同条件对照做一个审计面板，分开记录宿主溶血、细胞毒性、神经效应及物种。完整毒蛋白片段、前体和类似物单列，不作为已实测的短肽阳性。
2. 在同实验内查原六种子对连续测量的排序，并将同一骨架不同修饰的成对实验单列。排序失败与同输入下无法区分分子是两个不同问题。
3. 等面板中有足够的家族独立对照，再决定是否补同终点的 Cys=0 数据或加真实分子输入。旧溶血辅助实验在目标候选组漏判从 41 增至 43，支持先区分样本对象和终点。 [旧溶血辅助对照](../hemolysis_auxiliary_classifier_pilot_20261026/REPORT.md)

[67 条逐样本核查](逐样本核查.md) · [逐样本模型概率及证据 JSON](audit.json) · [复算脚本](build_audit.py)
