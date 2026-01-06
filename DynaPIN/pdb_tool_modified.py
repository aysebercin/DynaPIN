import os
import sys

def pad_line(line):
    """Helper function to pad line to 80 characters in case it is shorter"""
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    return line[:81]


def run(fhandle, chains, fh2):
    prev_line = None

    _pad_line = pad_line
    records = ('ATOM', 'HETATM', 'ANISOU')
    for line in fhandle:
        if line.startswith(records):
            line = _pad_line(line)
          
            segid = line[72:76].strip()[-1]
            if segid in chains:
                continue
            else:
                fh2.write(line[:21] + segid.ljust(1) + line[22:72] + '    ' + line[76:])
            prev_line = line
        elif line.startswith('TER'):  
            line = _pad_line(line)
            segid = prev_line[72:76].strip()[:-1]
            if segid in chains:
                continue
            else:
                fh2.write(line[:21] + segid.ljust(1) + line[22:72] + '    ' + line[76:])
        else:
            fh2.write(line)

def main2(file, chains):
    fh1 = open(file, 'r')
    fh2 = file.split('.')[0]+'_chained.pdb'
    with open(fh2, 'w+') as fhandle:
        run(fh1, chains, fhandle)

def main():
    pdbfh = open(sys.argv[1:][1], 'r')
    option = sys.argv[1:][0][1:]
    new_pdb = run(pdbfh, option)

    try:
        _buffer = []
        _buffer_size = 5000 
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                sys.stdout.write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        sys.stdout.write(''.join(_buffer))
        sys.stdout.flush()
    except IOError:
        pass

    pdbfh.close()
    sys.exit(0)
