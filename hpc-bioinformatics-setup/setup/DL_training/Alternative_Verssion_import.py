# In this case MRVIJAX from V 1.4 replaced by MRVI to use in V 1.3.3
import torch

model_path = '/work/archive/public_studies/merged_data/clean_dbl/mrvi_af_sex_model/model.pt'

checkpoint = torch.load(model_path, map_location='cpu', weights_only=False)

checkpoint['attr_dict']['registry_']['model_name'] = 'MRVI'
checkpoint['attr_dict']['registry_']['scvi_version'] = '1.3.3'

torch.save(checkpoint, model_path)

print(" Model metadata updated successfully!")

from scvi.external import MRVI
model = MRVI.load('/work/archive/public_studies/merged_data/clean_dbl/mrvi_af_sex_model', adata)
print(" Model loaded successfully!")
