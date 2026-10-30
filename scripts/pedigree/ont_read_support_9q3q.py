import subprocess, sys, json
S='/usr/bin/samtools'
D='/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/ont'
BAMS={'PAN011':[f'{D}/SRR3782605{i}.PAN011.v131.bam' for i in (4,5,6,7)],
      'PAN027':[f'{D}/SRR3782606{i}.PAN027.v131.bam' for i in (8,9)]}
LOCI=[
 ('PAN011','father','PAN011#2#chr9',133745035,'MODULE carrier (hap2)',133684561,133738660),
 ('PAN011','father','PAN011#1#chr9',136101078,'other hap (hap1)',None,None),
 ('PAN011','father','PAN011#2#chr3',205578940,'donor chr3 hap2',None,None),
 ('PAN011','father','PAN011#1#chr3',202554222,'other chr3 hap1',None,None),
 ('PAN027','SELF','PAN027#2#chr9',136209149,'transmitted pat chr9',136152460,136206559),
 ('PAN027','SELF','PAN027#1#chr9',129071589,'maternal chr9 (contrast)',None,None),
]
def stats(bams,contig,L,mod_lo,mod_hi,W=1000):
    def q(region):
        n=0
        for b in bams:
            p=subprocess.run([S,'view','-c',b,region],capture_output=True,text=True)
            try: n+=int(p.stdout.strip() or 0)
            except: pass
        return n
    term_lo=L-100_000
    term_n=q(f'{contig}:{term_lo}-{L}')
    tail_n=q(f'{contig}:{L-5_000}-{L}')
    out={'term_reads':term_n,'term_per_kb':round(term_n/100,2),'tail5kb_reads':tail_n}
    if mod_lo:
        mod_n=q(f'{contig}:{mod_lo}-{mod_hi}')
        out.update(module_reads=mod_n,module_per_kb=round(mod_n/((mod_hi-mod_lo)/1000),2))
    return out
res={}
for sample,who,contig,L,label,mlo,mhi in LOCI:
    r=stats(BAMS[sample],contig,L,mlo,mhi); r.update(who=who,label=label,contig=contig,length=L)
    res[f'{who}:{contig}']=r
    print(f"{who:5s} {contig:16s} {label:28s} {r}")
json.dump(res,open('read_support.json','w'),indent=1)
