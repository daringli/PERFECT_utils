def divide_interval_remainder(Npsi,N):
    """Divides an Npsi object into N piles.  If the piles cannot be evenly sized the last pile will have a different size, which will be smaller if possible."""
    remainder = Npsi%N
    if remainder == 0:
        # all will be equally sized
        return [Npsi//N]*N
    ret = [None] * N

    remainder = Npsi%(N-1)
    if remainder == 0:
        ret[:-1] = [Npsi//(N-1)-1]*(N-1)
        ret[-1] = N-1
    else:
        ret[:-1] = [Npsi//(N-1)]*(N-1)
        ret[-1] = remainder

    assert(sum(ret) == Npsi)

    return ret

if __name__ == "__main__":
    Npsi = 100
    a = range(0,Npsi)
    N = 5
    print divide_interval_remainder(Npsi,N)
    N = 6
    print divide_interval_remainder(Npsi,N)
    N = 7
    print divide_interval_remainder(Npsi,N)
    N = 8
    print divide_interval_remainder(Npsi,N)
    N = 9
    print divide_interval_remainder(Npsi,N)
    N = 10
    print divide_interval_remainder(Npsi,N)
    N = 11
    print divide_interval_remainder(Npsi,N)
    N = 12
    print divide_interval_remainder(Npsi,N)
    N = 13
    print divide_interval_remainder(Npsi,N)
