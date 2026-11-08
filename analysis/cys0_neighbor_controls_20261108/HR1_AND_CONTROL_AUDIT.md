# HR1 标签追溯和六条正确预测对照核查

## HR1

GHXTox 的 peptide_5088 裸序列为 INLKAIAALVKKVL，原标签 0。对应原 CAPTP 公开训练文件，六种子五折清单中 fold 1 和 fold 3 的 role=train，fold 0 的 role=test。旧工程已验证 CAPTP 公开训练文件的顺序、序列和标签与当前 6387 条一致。该公开文件没有逐条上游 accession 或标签实验来源。 [旧一致性记录](../ghxtox_error_atlas_20261025/evidence_phase2/captp_source_verification.json) · [CAPTP 数据](https://github.com/jiaoshihu/CAPTP/blob/main/data/train_data.txt)

2023 年原文表 1 中 HR1 的受试形式为 INLKAIAALVKKVL-NH2，人红细胞、320 µg/mL 溶血 56%，37 °C、1 小时。相同实验的 Mastoparan-L 为 100%。训练序列未标注末端，因此该实验支持对应酰胺分子的溶血，不能确定原标签依据的具体化学形式或终点。 [原文](https://www.mdpi.com/2072-6651/15/10/591)

上游精确匹配复用原有检索范围。CSM-Toxin 的两份公开阳性清单及对应 CD-HIT 清单、ToxinPred2 主阳性、ATSE 原训练和测试清单均未找到该裸序列。另核对 CSM-Toxin 仓库所含 ToxinPred2 主阴性文件，也未找到。CSM-Toxin UniProt 阴性原始大文件本轮请求未在时限内完成，CAPTP 单条来源仍未确认。不能把未命中已查清单解释为所有上游均不存在。

UniProt P17238 当前为带 toxin 等关键词的 60 aa 已审阅前体，包含该 14 aa 序列，后续有 Gly。该条说明天然分子的身份线索，不代表 2023 年 CAPTP 筛选时的记录状态或证明 HR1 来源于该前体。CAPTP 论文描述部分新增负类使用 UniProt 缺少毒素关键词的规则，实际单条是从哪个集合进入仍未知。 [UniProt P17238](https://www.uniprot.org/uniprotkb/P17238/entry) · [CAPTP 数据方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC11654579/)

## 六条正确预测对照

| 样本 | 旧对照分数 | 可核查的原始对象和实验 | 作为人红细胞对照的状态 |
|---|---:|---|---|
| peptide_357 | 0.944 | P86815 对应蛙皮分泌物缓激肽相关肽。UniProt 的平滑肌作用标记为相似性推断，本轮未定位该分子的直接实验原文。 | 待核原文 |
| peptide_290 | 0.935 | Cupiennin-6d，FINTIKLLIEKYREWKNKQSS-NH2。2012 年原文图 8 鉴定序列与 C 端酰胺形式。 | 该原文为毒液组成鉴定，未找到逐肽人红细胞实验 |
| peptide_1601 | 0.927 | Cupiennin-6c，FINTIRLLINKYREWKNKKDS-NH2。2012 年原文图 8 鉴定序列与 C 端酰胺形式。 | 同上 |
| peptide_280 | 0.955 | Cupiennin-5a，又称 SCP-7a，KFGKVLKFLAKTLAKHLAKKQAQS。2012 年原文图 8 鉴定该序列，图中未标 C 端酰胺。 | 同上 |
| peptide_720 | 0.924 | UniProt Q7M3T3 是 Vespa mandarinia 的 Vespakinin-M，记录第 4 位 4-羟脯氨酸。另一篇 Vespa magnifica 论文也称 Vespakinin-M，给出相同裸序列，但明确与含羟脯氨酸者有差异。 | 两个同名分子分别保留，本轮未找到当前 Q7M3T3 分子的同条件人红细胞实验 |
| peptide_774 | 0.964 | Vespakinin-X，ARPPGFSPFRIV。1977 年从 Vespa xanthoptera 毒液分离和测序，论文研究平滑肌、血压及毒液分离组分。 | 红细胞结果不能直接从粗毒液或其他组分归给该肽 |

蜘蛛肽来源为 [Trachsel 等 2012，图 8](https://doi.org/10.1111/j.1742-4658.2012.08650.x)。Vespakinin-X 来源为 [Yasuhara 等 1977](https://doi.org/10.1248/cpb.25.936)。Vespakinin-M 来源为 [1976 原论文书目](https://pubmed.ncbi.nlm.nih.gov/1017116/) 和 [Vespa magnifica 2006 原文](https://pubmed.ncbi.nlm.nih.gov/16870147/)。六条均尚缺与 HR1 或 Mastoparan-L 同物种、剂量和测量终点的直接溶血对照，未观察到实验不等于实验阴性。

## 下一步

优先在已核查的 Mastoparan 55 肽面板中寻找已有六种子折外预测且实验条件一致的肽，保留所有原标签而非只看正确预测的有毒类。配对样本应区分不同外折预测和模型，用同一外折模型评分的候选需事先核实。对 HR1 若能找到 CAPTP 数据生成时的逐条 accession、2023_04 UniProt 快照或上游阴性来源，可继续确定负标签从何而来。研究结论继续限定于溶血实验与原广义标签的关系。
