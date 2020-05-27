from glob import glob as local_glob
import os
from subprocess import check_output
import tempfile

def ssh_args():
  args = []
  if 'SSH_CRED' in os.environ:
    args += ['-i',os.environ['SSH_CRED']]
  return args

def _scp_cmd(src,dst):
  res = ['scp']+ssh_args()+[src,dst]
  #print("CMD:"+(' '.join(res)))
  return res

def is_ssh_path(fn):
  return ('@' in fn) and (':' in fn)

def local_ver(fn:str)->str:
  if is_ssh_path(fn):
    tmp_fn = tempfile.mktemp(prefix='wald',suffix='.'+fn.split('.')[-1])
    cmd = _scp_cmd(fn,tmp_fn)
    check_output(cmd)

    if fn.endswith('.hdf'):
      cmd = _scp_cmd(fn+'.xml',tmp_fn+'.xml')
      check_output(cmd)

    return tmp_fn

  return fn

def glob(pattern):
  if is_ssh_path(pattern):
    return remote_glob(pattern)
  return local_glob(pattern)

def remote_glob(pattern):
  conn,path = pattern.split(':')
  # print(conn,path)
  res = _remote_glob(conn,path)
  res = [r.rstrip('/') for r in res]
  return [f'{conn}:{p}' for p in res]

def _remote_glob(conn,path):
  components = path.split('/')
  static_path = ''
  wild_card=''
  for i,c in enumerate(components):
    if '*' in c:
      wild_card = c
      break
    static_path += c +'/'

  # print(static_path,wild_card)
  # print(i,len(components))
  cmd = ['ssh'] + ssh_args() + [conn,'ls','-d','--file-type',static_path+wild_card]
  # print('Running: %s'%(' '.join(cmd)))
  try:
    res = check_output(cmd).decode('utf-8').splitlines()
  except:
    res = []
  # print(res)
  res = [os.path.join(static_path,r) for r in res]
  if i < (len(components)-1):
    more_results = [_remote_glob(conn,os.path.join(r,*components[i+1:])) for r in res if r.endswith('/')]
    res = [r for mr in more_results for r in mr]
  # print(res[-10:])
  # print(len(res))
  return res


