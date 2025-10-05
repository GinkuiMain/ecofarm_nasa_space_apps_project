# Conversão mm/hr → mm/mês para diferentes comprimentos de mês
# (c) Você pode rodar como script Python puro.
# Obs.: usei Decimal para evitar erros de ponto flutuante.

from decimal import Decimal, getcontext
getcontext().prec = 50  # alta precisão

import pandas as pd

def mmhr_to_mmmes(value_mm_per_hr, days):
    """Converte mm/hr para mm/mês, dado o número de dias do mês."""
    hours = Decimal(24) * Decimal(days)
    return Decimal(value_mm_per_hr) * hours

def calc_all(label, v):
    """Gera todas as conversões para 28, 29, 30, 31 dias e mês médio (30.4375)."""
    results = {
        "label": label,
        "value_mm/hr": str(v),
        "28_days_mm/month": str(mmhr_to_mmmes(v, 28)),
        "29_days_mm/month": str(mmhr_to_mmmes(v, 29)),
        "30_days_mm/month": str(mmhr_to_mmmes(v, 30)),
        "31_days_mm/month": str(mmhr_to_mmmes(v, 31)),
        "mean_month_30.4375_days_mm/month": str(mmhr_to_mmmes(v, Decimal('30.4375'))),
    }
    return results

# --- Entradas ---

# Interpretação 1: string com pontos como separadores de milhar
raw_str = "900.176.540.017.128"
v1 = Decimal(raw_str.replace(".", ""))  # trata como 900176540017128 mm/hr (irreal, mas calculado a pedido)

# Interpretação 2: valor plausível lido do .nc (ex.: 0.090 mm/hr)
v2 = Decimal("0.090")

# --- Processamento ---

rows = []
rows.append(calc_all("Interpretation_1 (treated as 900,176,540,017,128 mm/hr)", v1))
rows.append(calc_all("Interpretation_2 (0.090 mm/hr)", v2))

df = pd.DataFrame(rows, columns=[
    "label", "value_mm/hr",
    "28_days_mm/month", "29_days_mm/month", "30_days_mm/month",
    "31_days_mm/month", "mean_month_30.4375_days_mm/month"
])

# Exibe resultado
print(df.to_string(index=False))

# (Opcional) salvar CSV
df.to_csv("conversao_mmhr_para_mmmes.csv", index=False)

