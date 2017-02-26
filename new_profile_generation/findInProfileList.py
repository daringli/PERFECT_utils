from profile import Profile

def findInProfileList(profileList,profile,species,generator="any",warnMultiple=True):
    if type(profileList) is not list:
        print "findInProfileList : error: profileList should be a list"
        raise ValueError('Profile list is not list')

    matches = []
    for i,p in enumerate(profileList):
        if type(p) is not Profile:
            print "findInProfileList : error: element <" + str(i) + "> in profile list not profile"
            #raise ValueError('Invalid element in profileList')
        
        #print p.species is species
        #print generator
        
        if (p.profile == profile or profile == "any") and (p.species == species or species == "any") and (p.generator == generator or generator == "any"):
            matches.append(i)

    if len(matches) == 0:
        print "findInProfileList : error: saught profile not in profileList (" + profile + ", " + species + ", " + generator + ")"
        raise ValueError('No matches in profileList')
    elif len(matches) == 1:
        return matches[0]
    elif warnMultiple:
        print "findInProfileList : error: more than one profile in the list of the type (" + profile + ", " + species + ", " + generator + ")"
        raise ValueError('Multiple matches in profileList')
    else:
        return matches

def extractProfilesFromList(profileList,wantedProfileList,wantedSpeciesList,wantedGeneratorList="any"):
    # verify inputs
    
    if type(wantedProfileList) is not list:
        wantedProfileList = [wantedProfileList]
    if type(wantedSpeciesList) is not list:
        wantedSpeciesList = [wantedSpeciesList]
    if type(wantedGeneratorList) is not list:
        wantedGeneratorList = [wantedGeneratorList]

    # if we want any species or profile, we are fine with multiple hits
    # but we need to remove other elements to avoid duplicates
    warnMultiple = True
    if any([e == "any" for e in wantedSpeciesList]):
        warnMultiple = False
        wantedSpeciesList = ["any"]
    if any([e == "any" for e in wantedProfileList]):
        warnMultiple = False
        wantedProfileList = ["any"]
    if any([e == "any" for e in wantedGeneratorList]):
        wantedGeneratorList = ["any"]

    
    indices = []
    for i in range(len(wantedProfileList)):
        if type(wantedProfileList[i]) is not str:
            print "findInProfileList : error: non-string element at position " + str(i) + " in wantedProfileList!"
            raise ValueError('Wanted profile not given as a string')
        for j in range(len(wantedSpeciesList)):
            if type(wantedSpeciesList[j]) is not str:
                print "findInProfileList : error: non-string element at position " + str(i) + " in wantedSpeciesList!"
                raise ValueError('Wanted species not given as a string')
            for k in range(len(wantedGeneratorList)):
                if type(wantedGeneratorList[k]) is not str:
                    print "findInProfileList : error: non-string element at position " + str(i) + " in wantedGeneratorList!"
                    raise ValueError('Wanted generator not given as a string')
                index = findInProfileList(profileList,wantedProfileList[i],wantedSpeciesList[j],wantedGeneratorList[k],warnMultiple=warnMultiple)
                if type(index) is list:
                    indices = indices + index
                else:
                    indices.append(index)
    return [profileList[index] for index in indices]

if __name__ == "__main__":
    class A(object):
        pass
    profileNames = ['n', 'ddx_n', 'eta', 'ddx_eta', 'T', 'ddx_T']
    specy = ''
    generator = ''
    Nprofiles = len(profileNames )
    profileList = [None]*Nprofiles
    for i in range(Nprofiles):
       profileList[i] = A()
       profileList[i].profile = profileNames[i]
       profileList[i].species = specy
       profileList[i].generator = generator
       
       
    print findInProfileList(profileList,"ddx_T",specy,generator="any") #3
    wantedProfiles = ["T","ddx_T"]
    [T,ddx_T] = extractProfilesFromList(profileList,wantedProfiles,specy)
    print T
    print ddx_T
