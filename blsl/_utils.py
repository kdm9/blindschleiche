import gzip

def rc(seq):
    d = {"A": "T", "G":"C", "C":"G", "T":"A"}
    return "".join([d.get(x, "N") for x in reversed(list(seq.upper()))])

def fqparse(stream, n=4):
    fqp = list()
    for line in stream:
        fqp.append(line.rstrip("\n"))
        if len(fqp) == n:
            yield fqp
            fqp = list()
    if len(fqp) == n:
        yield fqp
    assert len(fqp) == 0

