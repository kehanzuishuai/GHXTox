# 人红细胞末端修饰配对扩充

按同一基础序列、同一论文和实验批次、相同人红细胞及剂量、只改变一个末端状态的标准，从原论文核对溶血率。本轮新增 **3 篇原始研究、6 条序列、16 个剂量配对**。与既有记录合并后共有 **6 篇研究、10 条序列、23 个剂量配对**。其中 C 端游离酸／酰胺为 **4 篇、7 条序列、17 个配对**，其余为 N 端游离／乙酰化。浓度重复测量和同一序列的第二实验批次均保留为配对，样本量按序列另计。

| 肽 | 同批实验浓度 | C 端游离酸 | C 端酰胺 | 差值，百分点 |
| --- | ---: | ---: | ---: | ---: |
| Piscidin 3 | 100 µg/mL | 约 15% | 约 50% | +35 |
| Piscidin 1 | 100 µg/mL | 约 78% | 约 73% | −5，原文认为结果相近 |
| DMS-DA6 | 10 µM | 0.9% | 1.2% | +0.3 |
| DMS-DA6 | 50 µM | 7.4% | 3.7% | −3.7，报告波动较大 |
| PGLa | 32 µg/mL | 7% | 2% | −5 |
| PGLa | 1024 µg/mL | 10% | 61% | +51 |
| MAP | 32 µg/mL | 86% | 100% | +14 |

表中展示部分读数。PGLa 的负向读数保留，同时原文对整体剂量曲线的总结是游离酸形式溶血较低。DMS-DA6 的 50 µM 两组报告值分别为 7.4±5.5% 和 3.7±3.3%。这两种负差值用于记录实测趋势，现阶段缺少可证明稳定下降的独立复现。Piscidin 1 是作者明确描述为相近的对照。 [Piscidin 原文第 3.2 节](https://citeseerx.ist.psu.edu/document?doi=c27a0bbcfc094d6a3a58e09209a647764a829c5a&repid=rep1&type=pdf) · [Strandberg 原文表 3](https://publications.iupac.org/pac/pdf/2007/pdf/7904x0717.pdf) · [DMS-DA6 原文表 2](https://doi.org/10.1371/journal.pone.0205727.t002)

原文的浓度单位照录。两形式的 µg/mL 比较按原实验的名义质量浓度配对。MSI-103 两批分开，MAP 原表注明未实测的三个高浓度 100% 值剔除。旧资料里缺失的条件留空。台账保留报告的误差，避免把接近不变和实验不确定性当成等效证明。以直接溶血率作为定量主表，限值和效价进入附表。

| 额外筛查的研究 | 证据及处理 |
| --- | --- |
| AurH1 | 1000 µg/mL 两形式均低于 5%，组间差异未达到显著性。保留限值，文中 HC50 远超实测浓度。 [原文图 4](https://link.springer.com/article/10.1186/s12866-023-03090-7) |
| ChBac3.4、ChBac3.4-1 | 原表确认同序列 C 端配对，图 3 有低溶血曲线，待核逐点数据。溶血更高的 ChBac3.4-2 同时插入两个 Trp，不能归因于末端。 [原文](https://www.frontiersin.org/journals/cellular-and-infection-microbiology/articles/10.3389/fcimb.2020.552905/full) |
| Osmin | 人红细胞配对为 5.0／16.6 HU/mg，单列效价。 [原文第 3.5 节](https://www.academia.edu/22675548/Structural_identification_by_mass_spectrometry_of_a_novel_antimicrobial_peptide_from_the_venom_of_the_solitary_bee_Osmia_rufa_Hymenoptera_Megachilidae_) |
| Pteroicidin α | 图 4 有人红细胞曲线，待核逐点数据。旧库鱼红细胞值保留在历史记录。 [原论文收录于作者学位论文](https://sygal.normandie-univ.fr/informations/fichiers/telecharger/303) |
| Grammistin | 出版社、学会期刊目录及作者论文页均只找到摘要或获取全文入口，五条配对仍缺同条件数值。 [出版社](https://link.springer.com/article/10.1111/j.1444-2906.2006.01264.x) · [学会目录](https://www.miyagi.kopas.co.jp/JSFS/PUBS/FS/contents/72-5.html) |

与六种子训练清单按裸序列比对，**Piscidin 3 和 PGLa** 重合。它们可用于诊断序列分数为何无法区分两种分子形式。做新家族留出验证时，需把同研究、同家族和同骨架的样本成组考虑。当前 7 条 C 端序列适合诊断和设计对照，训练前优先补重复实验支持的下降样本和独立来源的低差异对照。

[合并后配对数据](human_rbc_matched_pairs.json) · [新增原始数据与出处](human_rbc_primary_sources.json) · [生成和检查脚本](build_human_rbc_panel.py)

下一轮补核见[不同家族的下降与低反应案例](HUMAN_RBC_FAMILY_EXTENSION.md)。
