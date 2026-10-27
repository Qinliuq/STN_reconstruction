import json,os,numpy as np,runpy,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
spec=runpy.run_path('sweeps/rescue_v2.py'); out='results/rescue_v2'
def load(L): return [json.load(open(f'sweeps/out/rescue_v2/{L}_s{s}/model_output_metrics.json')) for s in spec['seeds']]
def done(L): return all(os.path.exists(f'sweeps/out/rescue_v2/{L}_s{s}/model_output_metrics.json') for s in spec['seeds'])
def arr(ms):
    i=ms[0]['profile']['positions'].index(0.65)
    return {k:np.array(v) for k,v in dict(beta_pow=[m['beta_pow'] for m in ms],beta_frac=[m['beta_frac'] for m in ms],
      beta_frac65=[m['profile']['beta_frac'][i] for m in ms],sync=[m['sync_corr10ms'] for m in ms],sync65=[m['profile']['sync'][i] for m in ms],
      bursts=[m['frac_spikes_in_bursts'] for m in ms],pair_PVP=[m['beta_pow_pair_PVP_pop'] for m in ms],pair_PVN=[m['beta_pow_pair_PVN_pop'] for m in ms],
      sync_PVP=[m['sync_PVP_pop'] for m in ms],rate=[(12*m['rate_PVP_pop']+46*m['rate_PVN_pop'])/58 for m in ms],
      rate_PVP=[m['rate_PVP_pop'] for m in ms],rate_PVN=[m['rate_PVN_pop'] for m in ms],peak=[m['beta_peak_axis'] for m in ms]).items()}
h,u=arr(load('healthy')),arr(load('untreated'))
GAPS=['beta_pow','beta_frac','beta_frac65','sync','sync65','bursts','pair_PVP','pair_PVN','sync_PVP']
COMP=['beta_pow','beta_frac','sync65','bursts']
rows=[]
for c in spec['conditions']:
    if not done(c['label']): continue
    r=arr(load(c['label']))
    g={k:(u[k].mean()-r[k].mean())/(u[k].mean()-h[k].mean()) for k in GAPS}
    per=np.mean([(u[k]-r[k])/(u[k].mean()-h[k].mean()) for k in COMP],0)
    rows.append(dict(label=c['label'],group=c['group'],gaps=g,comp=per.mean(),comp_sem=per.std(ddof=1)/np.sqrt(len(per)),
                     dr=100*(r['rate'].mean()-u['rate'].mean())/u['rate'].mean(),rp=r['rate_PVP'].mean(),rn=r['rate_PVN'].mean(),
                     bu=r['bursts'].mean(),npk=int(np.sum((r['peak']>=0.5)&(r['peak']<=0.75)))))
with open(f'{out}/gap_table.csv','w') as fh:
    fh.write('label,'+','.join(f'gap_{k}' for k in GAPS)+',composite,composite_sem,rate_change_pct,rate_PVP,rate_PVN,bursts,beta_peak_upper_middle_seeds\n')
    for r in rows: fh.write(f"{r['label']},"+','.join(f"{r['gaps'][k]:.3f}" for k in GAPS)+f",{r['comp']:.3f},{r['comp_sem']:.3f},{r['dr']:.1f},{r['rp']:.2f},{r['rn']:.2f},{r['bu']:.4f},{r['npk']}\n")
# Kv3.1 dose-response
kv=[r for r in rows if r['group'] in ('untreated','Kv3.1')]
x=[1.0 if r['label']=='untreated' else float(r['label'][len('kv31x'):]) for r in kv]
fig,axs=plt.subplots(1,3,figsize=(15,4))
for k,lab in [('beta_pow','beta power'),('beta_frac','beta fraction'),('sync65','synchrony @0.65'),('bursts','bursts'),('pair_PVP','shared beta PV+ pairs'),('pair_PVN','shared beta PV- pairs')]:
    axs[0].plot(x,[r['gaps'][k] for r in kv],'o-',label=lab)
axs[0].axhline(1,color='k',ls='--',lw=0.8); axs[0].axhline(0,color='tab:red',ls='--',lw=0.8); axs[0].set_xlabel('Kv3.1 factor'); axs[0].set_title('gap closed per metric'); axs[0].legend(fontsize=7)
axs[1].errorbar(x,[r['comp'] for r in kv],yerr=[r['comp_sem'] for r in kv],fmt='o-',color='k',capsize=3)
axs[1].axhline(1,color='k',ls='--',lw=0.8); axs[1].axhline(0,color='tab:red',ls='--',lw=0.8); axs[1].set_xlabel('Kv3.1 factor')
axs[1].set_title('composite gap closed (beta pow, beta frac, sync@0.65, bursts)\nmean ± SEM over seeds',fontsize=9)
axs[2].plot(x,[r['rp'] for r in kv],'o-g',label='PV+'); axs[2].plot(x,[r['rn'] for r in kv],'o-b',label='PV-')
axs[2].axhline(h['rate_PVP'].mean(),color='g',ls=':',label='healthy PV+'); axs[2].axhline(h['rate_PVN'].mean(),color='b',ls=':',label='healthy PV-')
axs[2].set_xlabel('Kv3.1 factor'); axs[2].set_ylabel('rate (Hz)'); axs[2].set_title('rates'); axs[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{out}/kv31_dose.png',dpi=110)
# composite vs rate change
col={'CaT':'tab:orange','Kir':'tab:purple','HCN':'tab:green','Kv3.1':'tab:blue','with Kv3.1 x1.25':'tab:red','with Kv3.1 x1.3':'tab:brown'}
fig,axs=plt.subplots(1,3,figsize=(17,4.5))
for r in rows:
    if r['group'] not in col: continue
    for ax,v in zip(axs,[r['comp'],r['gaps']['bursts'],r['gaps']['beta_frac']]):
        ax.scatter(r['dr'],v,color=col[r['group']],s=35); ax.annotate(r['label'],(r['dr'],v),fontsize=6.5,xytext=(3,2),textcoords='offset points')
for ax,t in zip(axs,['composite gap closed','bursts gap closed','beta fraction gap closed']):
    ax.axhline(1,color='k',ls='--',lw=0.8); ax.axhline(0,color='tab:red',ls='--',lw=0.8); ax.set_xlabel('mean STN rate change vs untreated (%)'); ax.set_title(t)
handles=[plt.Line2D([],[],marker='o',ls='',color=c,label=g) for g,c in col.items()]; axs[0].legend(handles=handles,fontsize=7)
fig.tight_layout(); fig.savefig(f'{out}/gap_vs_rate.png',dpi=110)
for r in rows: print(f"{r['label']:18s} comp {r['comp']:5.2f}±{r['comp_sem']:.2f} dr {r['dr']:5.0f}")
