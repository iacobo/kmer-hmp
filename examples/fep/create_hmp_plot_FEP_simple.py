record_ids = sorted(list(set(Path(record.id).parent.stem for alignment in alignments for record in alignment)))
record_dict = {record_id:i for i, record_id in enumerate(record_ids)}
record_dict_reverse = {i:record_id for i, record_id in enumerate(record_ids)}

k = len(list(kmer_pvalues.keys())[0])

df_score = alignment2df(alignments,k,kmer_pvalues,record_dict)

num_tests = len(df_score)
thresh_amount = int(0.01*num_tests)
df_temp = df_score.nsmallest(thresh_amount, 'pval')

plt.scatter(df_temp.index.get_level_values('absolute_pos'), -np.log10(df_temp['pval']*len(df_score)))
