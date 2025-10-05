import pandas as pd

df = pd.read_csv('../df/SMAP_L3_SM_P_E_20251003_R19240_001.csv')

# verifica nulos
print(df.isnull().sum())

# Remove linhas com valores NaNs e gera novo CSV limpo
df = df.dropna()

# Mantém apenas uma amostra (10% da base original) para facilitar a plotagem
df = df.sample(frac=0.1, random_state=1)

df.to_csv('../df/SMAP_L3_SM_P_E_20251003_R19240_001_cleaned.csv', index=False)

# Analisa CSV limpo
df = pd.read_csv('../df/SMAP_L3_SM_P_E_20251003_R19240_001_cleaned.csv')
print(df.describe())

# retorna dados de Sergipe
LAT_MIN = -11.5
LAT_MAX = -9.5
LON_MIN = -38.5
LON_MAX = -36.0

sergipe_mask = (
    (df["Latitude"] >= LAT_MIN) &
    (df["Latitude"] <= LAT_MAX) &
    (df["Longitude"] >= LON_MIN) &
    (df["Longitude"] <= LON_MAX)
)

df_sergipe = df[sergipe_mask]

print(df_sergipe.describe())