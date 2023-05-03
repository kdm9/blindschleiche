
def rc(seq):
    d = {"A": "T", "G":"C", "C":"G", "T":"A"}
    return "".join([d.get(x, "N") for x in reversed(list(seq.upper()))])

def fqpair(stream):
    fqp = list()
    for line in stream:
        fqp.append(line.rstrip("\n"))
        if len(fqp) == 8:
            yield fqp
            fqp = list()
    if len(fqp) == 8:
        yield fqp
    else:
        assert len(fqp) == 0

def fqparse(stream):
    fqp = list()
    for line in stream:
        fqp.append(line.rstrip("\n"))
        if len(fqp) == 4:
            yield fqp
            fqp = list()
    if len(fqp) == 4:
        yield fqp
    else:
        assert len(fqp) == 0




