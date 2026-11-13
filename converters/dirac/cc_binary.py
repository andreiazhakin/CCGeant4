"""Portable, self-describing little-endian numeric arrays for CC tables v1."""
from pathlib import Path
import math
import os
import struct
import numpy as np

HEADER = struct.Struct('<8sIIIIQQ216s')
DESCRIPTOR = struct.Struct('<64s4sI6Q8s')
MAGIC = b'CCPACK01'
KINDS = {'delta': 1, 'global': 2, 'stopping': 3}
DTYPES = {'f8': np.dtype('<f8'), 'i8': np.dtype('<i8'), 'i4': np.dtype('<i4')}


def align(n):
    return (int(n)+63)//64*64


def write_pack(path, kind, zp, arrays):
    path = Path(path)
    if kind not in KINDS or not 1 <= len(arrays) <= 64:
        raise ValueError('Invalid table kind/array count')
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    prepared = []
    cursor = align(HEADER.size + len(arrays)*DESCRIPTOR.size)
    for name, value in arrays.items():
        if not name or len(name.encode('ascii')) >= 64:
            raise ValueError('Invalid array name')
        a = np.asarray(value)
        tag = a.dtype.kind+str(a.dtype.itemsize)
        if tag not in DTYPES or not 1 <= a.ndim <= 4 or not a.size:
            raise ValueError(f'Unsupported dtype/shape: {name} {a.dtype} {a.shape}')
        a = np.ascontiguousarray(a, dtype=DTYPES[tag])
        prepared.append((name, a, tag, cursor))
        cursor = align(cursor+a.nbytes)
    header_bytes = align(HEADER.size+len(prepared)*DESCRIPTOR.size)
    temporary = path.with_suffix(path.suffix+'.partial')
    with temporary.open('xb') as stream:
        stream.write(HEADER.pack(MAGIC, 1, KINDS[kind], int(zp), len(prepared),
                                 header_bytes, cursor, bytes(216)))
        for name, a, tag, offset in prepared:
            shape = list(a.shape)+[0]*(4-a.ndim)
            stream.write(DESCRIPTOR.pack(name.encode(), tag.encode(), a.ndim,
                                        *shape, offset, a.nbytes, bytes(8)))
        for name, a, tag, offset in prepared:
            stream.write(bytes(offset-stream.tell()))
            a.tofile(stream)
        stream.write(bytes(cursor-stream.tell()))
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)
    return {'path': str(path), 'bytes': path.stat().st_size,
            'kind': kind, 'Zp': int(zp),
            'arrays': {name: {'dtype': tag, 'shape': list(a.shape), 'offset': offset}
                       for name, a, tag, offset in prepared}}


def read_pack(path):
    path = Path(path)
    result = {}
    with path.open('rb') as stream:
        raw = stream.read(HEADER.size)
        if len(raw) != HEADER.size:
            raise ValueError('Truncated header')
        magic, version, kind, zp, count, header_bytes, size, reserved = HEADER.unpack(raw)
        if (magic != MAGIC or version != 1 or kind not in KINDS.values()
                or not 1 <= count <= 64 or size != path.stat().st_size
                or header_bytes != align(256+128*count) or any(reserved)):
            raise ValueError(f'Invalid CC table header: {path}')
        previous_end = header_bytes
        descriptors = []
        for _ in range(count):
            raw = stream.read(DESCRIPTOR.size)
            if len(raw) != DESCRIPTOR.size:
                raise ValueError('Truncated descriptor')
            name, tag, ndim, *rest = DESCRIPTOR.unpack(raw)
            name = name.rstrip(b'\0').decode('ascii')
            tag = tag.rstrip(b'\0').decode('ascii')
            shape, offset, nbytes, pad = rest[:4], rest[4], rest[5], rest[6]
            if (not name or name in result or tag not in DTYPES or not 1 <= ndim <= 4
                    or any(shape[ndim:]) or any(v <= 0 for v in shape[:ndim])
                    or offset % 64 or offset < previous_end or any(pad)
                    or nbytes != math.prod(shape[:ndim])*DTYPES[tag].itemsize
                    or offset+nbytes > size):
                raise ValueError(f'Invalid array descriptor: {path} {name}')
            result[name] = None
            descriptors.append((name, tag, tuple(shape[:ndim]), offset))
            previous_end = offset+nbytes
    for name, tag, shape, offset in descriptors:
        result[name] = np.memmap(path, mode='r', dtype=DTYPES[tag], offset=offset, shape=shape)
    return {'kind': kind, 'Zp': zp}, result
