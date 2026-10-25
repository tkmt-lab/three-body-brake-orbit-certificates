#!/usr/bin/env python3
"""Distribution integrity and bounded text/PDF privacy checks; no mathematical proof."""
import argparse,hashlib,json,re,subprocess
from verify import ROOT,load,require,A_ROWS,B_ROWS,C_ROWS,TWO

def files():
 return sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(x in ('.git','__pycache__','work') for x in p.relative_to(ROOT).parts) and p.name!='.DS_Store' and p.suffix!='.pyc' and p!=ROOT/'manifest.json')
def records():return [{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files()]
def run(write=False):
 actual=records()
 if write:
  (ROOT/'manifest.json').write_text(json.dumps({'version':'1.0.0','self_excluded':'manifest.json','excluded':['.git/','__pycache__/','*.pyc','work/','.DS_Store'],'files':actual},indent=2)+'\n');return {'status':'MANIFEST_WRITTEN','files':len(actual)}
 require(load('manifest.json')['files']==actual,'manifest membership/size/hash mismatch')
 for x in load('provenance.json')['files']:
  p=ROOT/x['path'];b=p.read_bytes();require(len(b)==x['bytes'] and hashlib.sha256(b).hexdigest()==x['sha256'],'provenance output binding '+x['path'])
 forbidden=[r'/'+s+r'/' for s in ['Users','Volumes','home','root']]
 forbidden+=[r'(?i)mail'+'box',r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}',r'(?i)PRIVATE'+' KEY',r'(?i)un'+'published',r'(?i)pre'+'publication',r'(?i)Version 1\.2',r'(?i)\(dra'+r'ft\)']
 pattern=re.compile('|'.join('(?:'+x+')' if not x.startswith('(?i)') else '(?i:'+x[4:]+')' for x in forbidden))
 for p in files():
  if p.suffix=='.pdf':continue
  text=p.read_text();require(not pattern.search(text),'privacy or superseded version text: '+str(p.relative_to(ROOT)))
 pdf=ROOT/'paper/article.pdf';text=subprocess.check_output(['pdftotext',str(pdf),'-'],text=True);info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
 require(not pattern.search(text+info),'PDF privacy/version');require('Version 1.0.0' in text,'PDF release version')
 require(re.search(r'Pages:\s+13',info),'expected 13-page reviewed layout')
 tex=(ROOT/'paper/article.tex').read_text();require(tex.count('VERIFICATION-RECORD.md')==1,'single supplement reference')
 require(len(re.findall(r'\\label\{([^}]+)\}',tex))==38,'paper labels')
 for token in ['1.0.0']:
  require(token in (ROOT/'CITATION.cff').read_text() and load('.zenodo.json')['version']==token,'metadata version')
 claims=load('claims.json');require(claims['same_zero_rows']==TWO and claims['certificate_count']==23 and claims['lean']['level']=='L0','scope metadata')
 require('unconfirmed' in claims['catalogue']['scope'] or 'nonmatches' in claims['catalogue']['scope'],'catalogue limitation')
 for family,rows in [('A',A_ROWS),('B',B_ROWS),('C',C_ROWS)]:require(len([p for p in (ROOT/'data/certificates').glob(f'{family}*.json')])==len(rows),'certificate membership')
 return {'status':'PASS','files':len(actual),'bytes':sum(x['bytes'] for x in actual),'pdf_pages':13,'privacy':'UTF-8 distribution text, PDF metadata and extracted text checked; bibliographic URLs are allowed','scope':'integrity and document metadata only'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write-manifest',action='store_true');a=p.parse_args();print(json.dumps(run(a.write_manifest),indent=2))
