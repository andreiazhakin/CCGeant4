"""Add the ATIMA source target molar mass to a five-column 2--40 MeV/u CSV.

Read masses from elemental spline headers and preserve the existing fields.
"""
from pathlib import Path
import argparse,gzip,json,math,re,tarfile

OLD='Zp,Zt,E_MeV_u,dEdx_MeV_cm2_g,dEdx_left_MeV_cm2_g'
NEW='Zp,Zt,At_ref_g_mol,E_MeV_u,dEdx_MeV_cm2_g,dEdx_left_MeV_cm2_g'

def convert(source,splines,output):
    if output.exists():raise FileExistsError(f'Refusing to overwrite {output}')
    pairs=set();rows=0
    with source.open(encoding='utf-8-sig') as stream:
        if stream.readline().strip()!=OLD:raise ValueError('Expected original five-column low-ATIMA CSV')
        for line in stream:
            fields=line.strip().split(',')
            if len(fields)!=5:raise ValueError('Malformed source row')
            z,zt=int(fields[0]),int(fields[1]);e,s,sl=map(float,fields[2:])
            if not (29<=z<=96 and 1<=zt<=96 and zt not in (85,86,87) and 2<=e<=40
                    and all(math.isfinite(v) and v>0 for v in (s,sl))):raise ValueError('Invalid source row')
            pairs.add((z,zt));rows+=1
    reference={};members={}
    with tarfile.open(splines,'r|') as archive:
        for member in archive:
            if not member.isfile() or not re.fullmatch(r'[A-Z][a-z]?\d+-[A-Z][a-z]?\.dat\.gz',Path(member.name).name):continue
            stream=archive.extractfile(member)
            with gzip.GzipFile(fileobj=stream) as unzipped:
                head=[unzipped.readline().decode('ascii').strip() for _ in range(4)]
            p=head[0].split();components=head[1].split();target=head[2].split()
            if len(p)!=2 or len(components)!=1 or components[0]!='1' or len(target)!=4:continue
            pair=(int(p[0]),int(target[0]))
            if pair not in pairs:continue
            if float(target[1])!=1:raise ValueError(f'Non-elemental source {member.name}')
            mass=float(target[3])
            if not math.isfinite(mass) or mass<=0:raise ValueError('Invalid source target mass')
            if pair in reference:raise ValueError(f'Duplicate source pair {pair}')
            reference[pair]=mass;members[pair]=member.name
    if set(reference)!=pairs:raise ValueError(f'Missing original spline headers: {sorted(pairs-set(reference))[:10]}')
    output.parent.mkdir(parents=True,exist_ok=True)
    temporary=output.with_suffix(output.suffix+'.partial')
    with source.open(encoding='utf-8-sig') as src,temporary.open('w',encoding='ascii',newline='') as dst:
        next(src);dst.write(NEW+'\n')
        for line in src:
            fields=line.strip().split(',');mass=reference[(int(fields[0]),int(fields[1]))]
            dst.write(','.join(fields[:2]+[format(mass,'.17g')]+fields[2:])+'\n')
    temporary.replace(output)
    with source.open(encoding='utf-8-sig') as src,output.open() as dst:
        next(src);next(dst)
        for a,b in zip(src,dst,strict=True):
            fields=b.rstrip('\r\n').split(',')
            if ','.join(fields[:2]+fields[3:])!=a.rstrip('\r\n'):raise AssertionError('Source fields changed')
    report={'source_csv':str(source.resolve()),
      'source_splines':str(splines.resolve()),
      'output':str(output.resolve()),'rows':rows,'pairs':len(pairs),
      'original_columns_exact_roundtrip':True,'mass_source':'original ATIMA elemental target header',
      'references':[{'Zp':p[0],'Zt':p[1],'At_ref_g_mol':reference[p],'member':members[p]} for p in sorted(pairs)]}
    output.with_suffix('.provenance.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['output','rows','pairs','original_columns_exact_roundtrip']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--atima-splines',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    convert(a.source,a.atima_splines,a.output)
