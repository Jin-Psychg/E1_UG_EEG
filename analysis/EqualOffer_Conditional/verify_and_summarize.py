"""Verify conditional outputs and summarize sensitivity; no fitting or recoding."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd

out=Path(sys.argv[1]).resolve()
root=Path(__file__).resolve().parents[2]
primary=root/'results/EqualOffer_Aligned_20260928_v2'
r=pd.read_csv(out/'conditional_omnibus.csv').fillna({'deleted_id':''})
e=pd.read_csv(out/'conditional_probabilities.csv').fillna({'deleted_id':''})
roster=pd.read_csv(out/'model_roster.csv').fillna({'deleted_id':''})
manifest=pd.read_csv(out/'input_manifest.csv')
for x in manifest.itertuples(): assert hashlib.sha256(Path(x.path).read_bytes()).hexdigest()==x.sha256,x.path
original=pd.read_csv(primary/'input_manifest.csv')
for x in original[original.role.str.endswith('_canonical_trials')].itertuples():
    assert hashlib.sha256(Path(x.path).read_bytes()).hexdigest()==x.sha256
assert np.allclose(np.exp(-r.Chisq/2)*(1+r.Chisq/2),r.p_raw,rtol=1e-10,atol=1e-14)
assert (r.df==4).all()
assert np.isfinite(e[['prob','asymp.LCL','asymp.UCL','logit_SE']]).all().all()
assert ((e['prob']>=0)&(e['prob']<=1)).all()
assert ((e['asymp.LCL']<=e['prob'])&(e['prob']<=e['asymp.UCL'])).all()
for label in ['primary','actor_selected','matched_no_actor']:
    a=r[r.model==label].sort_values('p_raw')
    expected=np.minimum.accumulate((a.p_raw.to_numpy()*2/np.arange(1,3))[::-1])[::-1].clip(max=1)
    assert np.allclose(expected,a.p_BH_two_experiments,rtol=1e-12)
counts=[]
for ex in ['E1','E2']:
    d=pd.read_csv(primary/f'{ex}_analysis_trials.csv')
    rawpath=original.loc[original.role==f'{ex}_canonical_trials','path'].iloc[0]
    raw=pd.read_csv(rawpath)
    kept=raw[raw.Offers_Other.isin([5,6])&raw.reaction.isin([1,2])&raw.RT.between(300,3000)]
    assert np.array_equal(kept[['participant_id','index']].values,d[['participant_id','index']].values)
    assert np.array_equal((kept.reaction==2).astype(int),d.reject_binary)
    des=pd.read_csv(primary/'descriptive_cells.csv')
    for em in ['neu','aff','dis','dom','enj']:
        z=d[(d.allocation=='5:5')&(d.emotion==em)]
        pp=z.groupby('participant_id').reject_binary.sum()
        a=des[(des.experiment==ex)&(des.allocation=='5:5')&(des.emotion==em)].iloc[0]
        assert len(z)==a.valid_trials and z.reject_binary.sum()==a.rejections
        counts.append(dict(experiment=ex,emotion=em,trials=len(z),rejections=int(z.reject_binary.sum()),rejectors=int((pp>0).sum()),largest_participant_count=int(pp.max()),largest_participant=str(pp.idxmax()),rejection_rate=float(z.reject_binary.mean())))
pd.DataFrame(counts).to_csv(out/'equal_offer_counts.csv',index=False)
summary=[]
for ex in ['E1','E2']:
    for label,base in [('primary_delete','primary'),('actor_delete','actor_selected')]:
        a=r[(r.experiment==ex)&(r.model==label)]
        planned=roster[(roster.experiment==ex)&(roster.model==label)]
        probs=e[(e.experiment==ex)&(e.model==label)]
        full=e[(e.experiment==ex)&(e.model==base)][['emotion','prob']]
        delta=probs.merge(full,on='emotion',suffixes=('','_full'))
        delta['absolute_pp_change']=100*abs(delta.prob-delta.prob_full)
        winner=delta.loc[delta.absolute_pp_change.idxmax()]
        summary.append(dict(experiment=ex,structure=base,accepted_deletions=len(a),planned_deletions=len(planned),missing_ids=';'.join(sorted(set(planned.deleted_id)-set(a.deleted_id))),p_min=float(a.p_raw.min()),p_max=float(a.p_raw.max()),largest_p_deleted_id=a.loc[a.p_raw.idxmax(),'deleted_id'],nominal_p_below_05=int((a.p_raw<.05).sum()),max_probability_change_pp=float(winner.absolute_pp_change),probability_change_deleted_id=winner.deleted_id,probability_change_emotion=winner.emotion))
pd.DataFrame(summary).to_csv(out/'deletion_summary.csv',index=False)
# Every literal reported numeric cell stays traceable to its output row and source model.
audit=[]
for file in ['conditional_omnibus.csv','conditional_probabilities.csv','equal_offer_counts.csv','deletion_summary.csv']:
    a=pd.read_csv(out/file,dtype=str,keep_default_na=False)
    for i,row in a.iterrows():
        for field,value in row.items():audit.append({'file':str(out/file),'row':i+2,'field':field,'literal':value})
pd.DataFrame(audit).to_csv(out/'stat_source_audit.csv',index=False)
receipt={'models_extracted':len(r),'source_hashes_rechecked':len(manifest),'canonical_trials_hash_and_rows_checked':True,'chisquare_p_independently_verified':True,'BH_two_test_families_verified':True,'probability_intervals_checked':True,'no_refits':True,'no_pairwise_significance_tests':True,'human_verified':False}
(out/'verification.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
(out/'STATUS.txt').write_text('COMPLETE: conditional arithmetic, input integrity and summaries verified.\n',encoding='utf-8')
print(r[r.model.isin(['primary','actor_selected','matched_no_actor'])][['experiment','model','Chisq','p_raw','p_BH_two_experiments']].to_string(index=False))
print(pd.DataFrame(summary).to_string(index=False))
print(e[e.model=='primary'][['experiment','emotion','prob','asymp.LCL','asymp.UCL']].to_string(index=False))

