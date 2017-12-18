def split_array(a,split_points):
    split_points.sort()

    if split_points[0] is not 0:
        split_points = [0] + split_points

    
    N = len(split_points) - 1
    b = [None]*N
    for i in range(0,N):
        if split_points[i+1] > len(a):
            split_points[i+1] = len(a)
        b[i] = a[split_points[i]:split_points[i+1]]
    return b

if __name__=="__main__":
    a = range(1,100)
    # normal use case
    split_points = [0,5,13,17,99]
    b = split_array(a,split_points)
    print b
    # unsorted
    split_points = [5,0,17,13,99]
    b = split_array(a,split_points)
    print b
    # unsorted, no leading zero after sort
    split_points = [5,3,13,17,99]
    b = split_array(a,split_points)
    print b
    # value larger than max
    split_points = [5,3,13,17,1010]
    b = split_array(a,split_points)
    print b
    # duplicate values after max has been floored
    split_points = [5,3,3,13,17,1010,1300]
    b = split_array(a,split_points)
    print b
