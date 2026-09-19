"""Verified, bounded-memory and resumable public-source downloads."""
import hashlib
import json
import os
import shutil
import ssl
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.request import Request, urlopen

CHUNK = 16 * 1024 * 1024

def verified_urlopen(request, timeout=90):
    """Use native certificate-chain validation without relaxing TLS verification."""
    import truststore
    context = truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    if not context.check_hostname or context.verify_mode != ssl.CERT_REQUIRED:
        raise RuntimeError("HTTPS certificate and hostname verification must remain enabled")
    return urlopen(request, timeout=timeout, context=context)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def dump(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    os.replace(tmp, path)

def download(spec, folder, workers=8):
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    path = folder / spec['file']
    size, digest = spec['bytes'], spec['sha256']
    if path.exists():
        if path.stat().st_size == size and sha(path) == digest:
            return path
        raise RuntimeError(f'Existing input failed checksum: {path}. Preserve it and report the error.')
    parts = folder / (path.name + '.parts'); parts.mkdir(exist_ok=True)
    ranges = [(i, min(i + CHUNK, size)) for i in range(0, size, CHUNK)]
    def valid(a, b):
        p = parts / str(a); marker = parts / (str(a) + '.sha256')
        return p.exists() and p.stat().st_size == b-a and marker.exists() and sha(p) == marker.read_text().strip()
    missing = [(a, b) for a, b in ranges if not valid(a, b)]
    required = size + sum(b-a for a,b in missing) + 1024**3
    if shutil.disk_usage(folder).free < required:
        raise RuntimeError(f'Insufficient free disk for {path.name}; need {required/1024**3:.1f} GiB.')
    if missing:
        # Normal public HTTP redirects only. A refusal is propagated, never bypassed.
        req = Request(spec['url'], headers={'Range':'bytes=0-0', 'Accept-Encoding':'identity'})
        with verified_urlopen(req, timeout=90) as r:
            resolved = r.url
        def part(ab):
            a,b=ab; p=parts/str(a)
            req=Request(resolved, headers={'Range':f'bytes={a}-{b-1}', 'Accept-Encoding':'identity'})
            with verified_urlopen(req, timeout=180) as r:
                expected=f'bytes {a}-{b-1}/{size}'
                if r.status != 206 or r.headers.get('Content-Range') != expected:
                    raise RuntimeError(f'Range not honored for {path.name}: {r.status}, {r.headers.get("Content-Range")}')
                tmp=p.with_suffix('.tmp'); h=hashlib.sha256(); n=0
                with tmp.open('wb') as f:
                    for data in iter(lambda:r.read(1024*1024), b''):
                        f.write(data); h.update(data); n+=len(data)
                if n != b-a: raise RuntimeError(f'Truncated range {a}: {n} != {b-a}')
            os.replace(tmp,p); p.with_name(p.name+'.sha256').write_text(h.hexdigest())
            return b-a
        pool=ThreadPoolExecutor(max_workers=workers)
        futures=[pool.submit(part,ab) for ab in missing]
        try:
            for i,future in enumerate(as_completed(futures),1):
                future.result()
                if i%16==0 or i==len(missing): print(f'{path.name}: {i}/{len(missing)} remaining chunks',flush=True)
        except BaseException as error:
            print(f'Download stopped: {error}. Cancelling pending chunks.',flush=True)
            for future in futures:future.cancel()
            raise
        finally:
            pool.shutdown(wait=True,cancel_futures=True)
    tmp=path.with_name(path.name+'.assembling')
    with tmp.open('wb') as out:
        for a,b in ranges:
            if not valid(a,b): raise RuntimeError('Chunk changed during assembly')
            with (parts/str(a)).open('rb') as src: shutil.copyfileobj(src,out,1024*1024)
    if tmp.stat().st_size != size or sha(tmp) != digest:
        raise RuntimeError(f'Full source checksum failed: {tmp}')
    os.replace(tmp,path); shutil.rmtree(parts)
    return path
