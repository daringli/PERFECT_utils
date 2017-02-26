def findUnsharedElements(list1,list2,list2GreaterThanList1Warning=False):
    """Returns a list of elements not found in both list1 and list2.
    Computes results by removing duplicate elements and makes no assumptions about the list being sorted.
    Example: list1 = [6,1,1,2,5,6] and list2 = [1,4,2,6,6] returns [1,5,4]"""

    # to not change the input lists
    list1=list(list1)
    list2=list(list2)
    
    if len(list1) < len(list2):
        list1,list2 = list2,list1
        
    i = 0
    while i< len(list1):
        delete = False
        j = 0
        while j < len(list2) and i< len(list1):
            if list1[i] == list2[j]:
                del list1[i]
                del list2[j]
                delete = True
            else:
                # only add to counter if we did not delete,
                # since list shrinks during deletion
                j=j + 1
        if not delete:
            i = i + 1

    return list1 + list2

if __name__ == "__main__":
    list1 = [6,1,1,2,5,6]
    list2 = [1,4,2,6,6]
    print findUnsharedElements(list1,list2) # [1,5,4]
    
    list1 = [1]
    list2 = [1]
    print findUnsharedElements(list1,list2) # []
    
    list1 = [1]
    list2 = []
    print findUnsharedElements(list1,list2) # [1]

    list1 = [6,4,5]
    list2 = [5,4,6,7]
    print findUnsharedElements(list1,list2) # [7]
    # we do not modify mutable inputs since we create new lists in func
    print list1
    print list2

    list1 = ['c','a']
    list2 = ['a','b']
    print findUnsharedElements(list1,list2) # ['c','b']
