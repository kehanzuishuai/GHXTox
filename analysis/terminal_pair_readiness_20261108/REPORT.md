# GHXTox 末端化学配对盘点

**这批数据适合做诊断。现有可核对数据不足以训练并独立验证“序列＋末端化学”的预测头。** 本轮沿用六种子 ProtT5 的原始序列评分作为历史参考，实测终点限定为红细胞溶血。

## 去重和证据分级

| 证据层级 | 裸序列数 | 可用范围 |
| --- | ---: | --- |
| 原表可核对的同条件定量数据 | 6 | 5 篇研究，19 个配对剂量点。人红细胞 4 条序列、3 篇研究，其中 C 端配对只有 Mast 21 |
| 原图支持的配对 | 1 | Pteroicidin α。人、海鲈和小斑猫鲨分别记录，旧库的通用 fish 数值对暂不纳入定量 |
| 摘要支持的配对 | 5 | Grammistin Gs1、Gs2、Pp1、Pp2b、Pp3，来自同一篇研究。待补原表和实验条件 |

共 12 条不同裸序列，包含旧 7 条和 2006 年 Grammistin 实验的 5 条。Gs2/Pp4a、Pp1/GsD、Pp2b/GsE 按别名去重。MreB1–9 与其 N 端加 W 的版本保留两条序列，验证时归同一设计组。Grammistin 五条来自同一研究，留研究验证时应一起留出。重复剂量属于同一肽骨架。19 个定量剂量点来自旧 19 个剂量块移出 Pteroicidin α 的 1 个，再加入 RR4 原表的 1 个。 [Grammistin 配对原文](https://link.springer.com/article/10.1111/j.1444-2906.2006.01264.x) · [序列原研究](https://www.jstage.jst.go.jp/article/fishsci1994/67/1/67_1_163/_pdf) · [别名研究](https://pubmed.ncbi.nlm.nih.gov/15777955/)

## 逐项核对

| 配对 | 比较 | 原始证据和处理 |
| --- | --- | --- |
| BP76 / BP77 | N 端游离／乙酰化 | 人红细胞，37 °C、1 h。表 1 的 34%／40% 对应 **150 µM**，旧库条目 1843、1847 写成 100 µM。本台账保留旧值和校正值。 [原表](https://www.academia.edu/18839892/A_library_of_linear_undecapeptides_with_bactericidal_activity_against_phytopathogenic_bacteria) |
| Mast 21 / Mast 21N | C 端游离酸／酰胺 | 人红细胞，37 °C、1 h、10 µM。溶血率 **5.4%／12.3%**。表 2 的脂质体泄漏率另列。 [原文表 1、2](https://www.jstage.jst.go.jp/article/bbb/64/5/64_5_985/_pdf) |
| MreB1–9 与 W-MreB1–9 | 分别比较 N 端游离／乙酰化 | 人红细胞，37 °C、1 h、50 和 100 µM。MreB1–9 两个浓度均报告 0%／0%，作为无可测差异的对照。W-MreB1–9 两个剂量都有分离。 [原文表 1、5](https://pmc.ncbi.nlm.nih.gov/articles/PMC5322399/) |
| Pteroicidin α | C 端游离酸／酰胺 | 图 4 支持三种红细胞分别配对。正文的海鲈游离酸 50% 与酰胺化后最高约 90% 不能直接当作同一物种同一剂量的精确数据点。 [原文书目](https://pubmed.ncbi.nlm.nih.gov/29108968/) · [作者学位论文收录原文，3.4 节及图 4](https://sygal.normandie-univ.fr/informations/fichiers/telecharger/303) |
| Kn2-5R | C 端游离酸／酰胺 | 羊红细胞，1–256 µM 的九个剂量分别配对。 [原文表 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC8348142/#molecules-26-04654-t004) |
| RR4 | C 端游离酸／酰胺 | 复核为马红细胞。旧 150、375 µM 配对外，补入原表 250 µM 的 9.8%／85.0%。 [原文及补充表 1、4](https://doi.org/10.3389/fmicb.2021.753874) |
| Grammistin 五条 | C 端游离酸／酰胺 | 摘要报告 Gs1 两个形式未观察到活性差异，其余四条酰胺化形式活性更高。待核红细胞来源、剂量、终点计量和分形式读数。 [原文摘要](https://link.springer.com/article/10.1111/j.1444-2906.2006.01264.x) |

## 对六种子困难样本的含义

- Gs2 和 Pp3 都与既有 Cys=0 共同漏判名单精确匹配，六种子 OOF 毒性分分别为 **0.0310** 和 **0.0357**。Gs1 也在原困难样本名单中，原论文报告其两种形式未观察到活性差异。这个家族内部的差异是选择对照的依据。Gs1、Gs2、Pp2b、Pp3 的裸序列都与六种子训练集清单重合，应作为历史诊断。 [逐样本历史核查](../cys0_common_error_assay_20261108/逐样本核查.md) · [Grammistin 原文摘要](https://link.springer.com/article/10.1111/j.1444-2906.2006.01264.x)
- 旧 19 个剂量块的每个配对在数据库中的原始 SMILES 都相同。台账用裸序列及明确区分的末端字段保留两种形式，Grammistin 摘要未确认的 N 端字段留空。相同 FASTA 对六种子序列模型必然给出相同的配对分数，说明输入存在信息缺口。现有数据还无法确定末端信息能救回多少共同漏判，也无法确定可靠结构有多少额外收益。 [历史原始字段核对](../hemolysis_three_direction_review_20261026/chemistry_metadata_parent_verification.json)

## 下一步

先补 Grammistin 2006 年全文，优先提取 Gs2、Pp3 及无明显差异对照 Gs1。并行补不同来源的人红细胞 C 端实测配对，保留无变化和反方向的对照。记录红细胞物种、剂量、实验时间、缓冲液和读数。固定剂量溶血率、HC50、其他效价和细胞活性分别保存。

有足够不同来源和家族的同终点数据后，先做按来源和骨架分组的留出比较，用同一批样本比较“序列”和“序列＋明确末端状态”的配对排序、溶血差值误差和误报。再在有可靠结构的同一批分子上测试结构的额外贡献。本轮只完成证据盘点，未训练模型或使用 Test1/Test2。

[配对台账](pair_ledger.json) · [生成脚本](build_ledger.py) · [Mast 21 原文](https://www.jstage.jst.go.jp/article/bbb/64/5/64_5_985/_pdf)

后续扩展见[跨论文人红细胞配对数据](HUMAN_RBC_PANEL.md)。
