def generic_species_labels(species_list,translate_dict={"D":"i", "He":"z","e":"e","N":"z"}):
    #replaces species names with other labels according to the supplied dictionary.
    for key in translate_dict.keys():
        if key in species_list:
            species_list[species_list.index(key)]=translate_dict[key]
    return species_list

if __name__=="__main__":
    sl=["D","e","He"]
    print generic_species_labels(sl)
    sl=["D","e","He"]
    print generic_species_labels(sl,{"He":"i"})
