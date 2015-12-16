import re,sys

def pretty_parse(dir):
    dir=dir.rstrip('/')
    line=dir.rsplit('/')[-1]
    print line
    match=re.match(r'Z([0-9]+)', line)
    if match:
        return r"$Z="+match.group(1)+"$"
    else:
        match=re.match(r'ZeffM1_([0-9]+\.?[0-9]*)', line)
        if match:
            print match.group(1)
            return r'$Z_{eff}-1='+match.group(1)+'$'
        else:
            match=re.match(r'n([0-9]+)', line)
            if match:
                if match.group(1)[0]=="0":
                    return "$n_d=0."+match.group(1)[1:]+"n_{d0}$"
                elif match.group(1)=="1":
                    return "$n=n_{d0}$"
                else:
                    return "$n_d="+match.group(1)+"n_{d0}$"
            else:
                match=re.match(r'T([0-9]+)', line)
                if match:
                    if match.group(1)[0]=="0":
                        return "$T=0."+match.group(1)[1:]+"T_{0}$"
                    elif match.group(1)=="1":
                        return "$T=T_{0}$"
                    else:
                        return "$T="+match.group(1)+"T_{0}$"
                else:
                    match=re.match(r'sp([A-Z][a-z]?)', line)
                    if match:
                        return match.group(1)
                    else:
                        match=re.match(r'([0-9]+\.?[0-9]*)', line)
                        if match:
                            return "$dn_d/dr="+match.group(1) +"\,m/s$"
                       

if __name__=='__main__':
    dirs=sys.argv[1:]
    print dirs
    pretty=[pretty_parse(dir) for dir in dirs]
    print pretty
