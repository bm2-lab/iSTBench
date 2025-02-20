## DR.S
file (abbreviation: f): The file path that stores the result of data consolidation when the number of domains is set. Refer to the example data for the specific format. This is a required parameter.
cluster_number (abbreviation: n): The number of classes of samples. This is a required parameter.
abundance_data (abbreviation: a): The file name of slice representation based on domain abundance. It should be stored in "file". This is a required parameter.
predicted_domain (abbreviation: p): The name of metadata that store identified domians. This is an optional parameter, with the default set to "predicted_domain".
slices_class (abbreviation: s): The name of metadata that store sample class. This is an optional parameter, with the default set to "slices_class".
sd_filter (abbreviation: S): The threshold for filtering slice representation based on domain-domain similarity. If the similarity between two domains is below the sd_filter value, their similarity will not be considered when representing the slice. This is an optional parameter, with the default set to 0.
