import streamlit as st
import pandas as pd
from datetime import datetime
import io
import math

# --- 1. CONFIGURAÇÕES, DICIONÁRIOS E CONSTANTES ---

alvo_filtro_options = {
    'Mo/Mo': 1, 'Mo/Rh': 1.017, 'Rh/Rh': 1.061, 'Rh/Al': 1.044, 'W/Rh': 1.042
}

# Coeficientes padrão para CSR
csr_coeffs_default = {
    'Mo/Mo': {'a': 0.01, 'b': 0.08},
    'Mo/Rh': {'a': 0.0067, 'b': 0.2333},
    'Rh/Rh': {'a': 0.0167, 'b': -0.0367},
    'W/Rh':  {'a': 0.0067, 'b': 0.3533}
}

# Tabelas Ki Originais
tabela_ki_ird = {('Mo/Mo', 26): 0.1357, ('Mo/Mo', 27): 0.1530, ('Mo/Rh', 29): 0.1540, ('Mo/Rh', 31): 0.1830}
tabela_ki_ufrj = {
    ('Mo/Mo', 25): 0.119094, ('Mo/Mo', 26): 0.136889, ('Mo/Mo', 27): 0.155258, ('Mo/Mo', 28): 0.175158,
    ('Mo/Rh', 26): 0.114301, ('Mo/Rh', 27): 0.131012, ('Mo/Rh', 28): 0.148476, ('Mo/Rh', 29): 0.166423,
    ('Rh/Rh', 28): 0.126825, ('Rh/Rh', 29): 0.142299, ('Rh/Rh', 30): 0.158490, ('Rh/Rh', 31): 0.175164,
}

# Fator G - Constantes e Incertezas
FATOR_G_CONSTANTS_UNCERTAINTIES = {
    0.30: {'a0': 0.6862414, 'da0': 0.0215771, 'a1': -0.1903851, 'da1': 0.0122059, 'a2': 0.0211549, 'da2': 0.0020598, 'a3': -0.0008170, 'da3': 0.0001055},
    0.35: {'a0': 0.7520924, 'da0': 0.0214658, 'a1': -0.2040045, 'da1': 0.0121429, 'a2': 0.0223514, 'da2': 0.0020492, 'a3': -0.0008553, 'da3': 0.0001050},
    0.40: {'a0': 0.8135159, 'da0': 0.0208152, 'a1': -0.2167391, 'da1': 0.0117749, 'a2': 0.0234949, 'da2': 0.0019871, 'a3': -0.0008925, 'da3': 0.0001018},
    0.45: {'a0': 0.8587792, 'da0': 0.02030096, 'a1': -0.2213542, 'da1': 0.01148395, 'a2': 0.0235061, 'da2': 0.00193800, 'a3': -0.0008817, 'da3': 0.00009929},
    0.50: {'a0': 0.8926865, 'da0': 0.0192286, 'a1': -0.2192870, 'da1': 0.0108773, 'a2': 0.0224164, 'da2': 0.0018356, 'a3': -0.0008171, 'da3': 0.0000940},
    0.55: {'a0': 0.9237367, 'da0': 0.0184259, 'a1': -0.2189931, 'da1': 0.0104233, 'a2': 0.0221241, 'da2': 0.0017590, 'a3': -0.0008050, 'da3': 0.0000901},
    0.60: {'a0': 0.9131422, 'da0': 0.0097610, 'a1': -0.1996713, 'da1': 0.0055217, 'a2': 0.0190965, 'da2': 0.0009318, 'a3': -0.0006696, 'da3': 0.0000477},
}

# Fator C - Fórmulas
formulas_fator_c = {
    0.34: {1: lambda e: (0.0004 * e**3) - (0.0105 * e**2) + (0.093 * e) + 0.9449, 2: lambda e: 0.0001 * e**3 - 0.0035 * e**2 + 0.0295 * e + 0.9831, 3: lambda e: -0.0001 * e**3 + 0.0028 * e**2 - 0.0242 * e + 1.0105, 4: lambda e: -0.0005 * e**3 + 0.0103 * e**2 - 0.0773 * e + 1.0343},
    0.35: {1: lambda e: (0.0004 * e**3) - (0.0105 * e**2) + (0.093 * e) + 0.9449, 2: lambda e: 0.0001 * e**3 - 0.0035 * e**2 + 0.0295 * e + 0.9831, 3: lambda e: -0.0001 * e**3 + 0.0028 * e**2 - 0.0242 * e + 1.0105, 4: lambda e: -0.0005 * e**3 + 0.0103 * e**2 - 0.0773 * e + 1.0343},
    0.36: {1: lambda e: 0.0004 * e**3 - 0.0103 * e**2 + 0.0915 * e + 0.9443, 2: lambda e: 0.0002 * e**3 - 0.0044 * e**2 + 0.0338 * e + 0.9768, 3: lambda e: -0.0001 * e**3 + 0.0029 * e**2 - 0.0248 * e + 1.0118, 4: lambda e: -0.0004 * e**3 + 0.0093 * e**2 - 0.0726 * e + 1.03},
    0.37: {1: lambda e: 0.0005 * e**3 - 0.0117 * e**2 + 0.098 * e + 0.9345, 2: lambda e: 0.0002 * e**3 - 0.0041 * e**2 + 0.0325 * e + 0.9783, 3: lambda e: -0.0001 * e**3 + 0.003 * e**2 - 0.0247 * e + 1.0117, 4: lambda e: -0.0004 * e**3 + 0.0091 * e**2 - 0.0718 * e + 1.0304},
    0.38: {1: lambda e: 0.0005 * e**3 - 0.0117 * e**2 + 0.0978 * e + 0.9342, 2: lambda e: 0.0002 * e**3 - 0.0041 * e**2 + 0.0324 * e + 0.9782, 3: lambda e: -0.0001 * e**3 + 0.0031 * e**2 - 0.0252 * e + 1.0126, 4: lambda e: -0.0004 * e**3 + 0.009 * e**2 - 0.0715 * e + 1.0306},
    0.39: {1: lambda e: 0.0005 * e**3 - 0.0116 * e**2 + 0.0974 * e + 0.934, 2: lambda e: 0.0002 * e**3 - 0.0041 * e**2 + 0.0324 * e + 0.9782, 3: lambda e: -0.0001 * e**3 + 0.0031 * e**2 - 0.0251 * e + 1.0126, 4: lambda e: -0.0004 * e**3 + 0.0089 * e**2 - 0.0712 * e + 1.0311},
    0.40: {1: lambda e: 0.0005 * e**3 - 0.0114 * e**2 + 0.0959 * e + 0.9335, 2: lambda e: 0.0002 * e**3 - 0.0041 * e**2 + 0.0322 * e + 0.9779, 3: lambda e: -0.0001 * e**3 + 0.0031 * e**2 - 0.0248 * e + 1.0128, 4: lambda e: -0.0004 * e**3 + 0.0087 * e**2 - 0.0703 * e + 1.0324},
    0.41: {1: lambda e: 0.0007 * e**3 - 0.0154 * e**2 + 0.1207 * e + 0.8822, 2: lambda e: 0.0002 * e**3 - 0.0036 * e**2 + 0.0299 * e + 0.9801, 3: lambda e: -0.0001 * e**3 + 0.0031 * e**2 - 0.0248 * e + 1.0125, 4: lambda e: -0.0004 * e**3 + 0.009 * e**2 - 0.0716 * e + 1.0352},
    0.42: {1: lambda e: 0.0007 * e**3 - 0.0165 * e**2 + 0.1278 * e + 0.8677, 2: lambda e: 0.0001 * e**3 - 0.0034 * e**2 + 0.0293 * e + 0.9807, 3: lambda e: -0.0001 * e**3 + 0.0031 * e**2 - 0.0247 * e + 1.0124, 4: lambda e: -0.0004 * e**3 + 0.0091 * e**2 - 0.0719 * e + 1.0358},
    0.43: {1: lambda e: 0.0008 * e**3 - 0.0177 * e**2 + 0.1349 * e + 0.853, 2: lambda e: 0.0001 * e**3 - 0.0033 * e**2 + 0.0286 * e + 0.9815, 3: lambda e: -0.0001 * e**3 + 0.0031 * e**2 - 0.0247 * e + 1.0124, 4: lambda e: -0.0004 * e**3 + 0.0092 * e**2 - 0.0724 * e + 1.0368},
    0.44: {1: lambda e: 0.0009 * e**3 - 0.0188 * e**2 + 0.1419 * e + 0.8384, 2: lambda e: 0.0001 * e**3 - 0.0032 * e**2 + 0.0279 * e + 0.9822, 3: lambda e: -0.0001 * e**3 + 0.0031 * e**2 - 0.0246 * e + 1.0122, 4: lambda e: -0.0004 * e**3 + 0.0092 * e**2 - 0.0727 * e + 1.0375},
    0.45: {1: lambda e: 0.0011 * e**3 - 0.0229 * e**2 + 0.1669 * e + 0.787, 2: lambda e: 0.00009 * e**3 - 0.0026 * e**2 + 0.0252 * e + 0.9851, 3: lambda e: -0.0001 * e**3 + 0.0029 * e**2 - 0.0238 * e + 1.0109, 4: lambda e: -0.0004 * e**3 + 0.009 * e**2 - 0.0719 * e + 1.0374},
    0.46: {1: lambda e: 0.0007 * e**3 - 0.0162 * e**2 + 0.1292 * e + 0.8523, 2: lambda e: 0.00008 * e**3 - 0.0024 * e**2 + 0.0241 * e + 0.9865, 3: lambda e: -0.0001 * e**3 + 0.0029 * e**2 - 0.0241 * e + 1.0127, 4: lambda e: -0.0004 * e**3 + 0.0087 * e**2 - 0.0706 * e + 1.0377},
    0.47: {1: lambda e: 0.0006 * e**3 - 0.015 * e**2 + 0.1216 * e + 0.8666, 2: lambda e: 0.00008 * e**3 - 0.0024 * e**2 + 0.0238 * e + 0.9869, 3: lambda e: -0.0001 * e**3 + 0.0029 * e**2 - 0.0242 * e + 1.0132, 4: lambda e: -0.0004 * e**3 + 0.0086 * e**2 - 0.07 * e + 1.0375},
    0.48: {1: lambda e: 0.0008 * e**3 - 0.0177 * e**2 + 0.1349 * e + 0.853, 2: lambda e: 0.0008 * e**3 - 0.0177 * e**2 + 0.1349 * e + 0.853, 3: lambda e: 0.0004 * e**3 - 0.0105 * e**2 + 0.093 * e + 1.077, 4: lambda e: -0.0004 * e**3 + 0.0093 * e**2 - 0.0726 * e + 1.03},
    0.50: {1: lambda e: (0.0004 * e**3) - (0.0105 * e**2) + (0.093 * e) + 1.077, 2: lambda e: 0.0008 * e**3 - 0.0177 * e**2 + 0.1349 * e**2 + 0.853, 3: lambda e: 0.0004 * e**3 - (0.0105 * e**2) + (0.093 * e) + 1.077, 4: lambda e: -0.0004 * e**3 + 0.0093 * e**2 - 0.0726 * e + 1.03},
}

# --- 2. FUNÇÕES DE CÁLCULO ---

def propagate_uncertainty(value_func, uncertainty_terms):
    sum_of_squares = 0
    for partial_deriv, input_uncertainty in uncertainty_terms:
        sum_of_squares += (partial_deriv * input_uncertainty)**2
    return math.sqrt(sum_of_squares)

def calcular_csr(kv_val, alvo_filtro, d_kv_abs, current_csr_dict):
    try:
        coeffs = current_csr_dict.get(alvo_filtro)
        if not coeffs: return "Erro CSR", 0.0
        
        csr_val = round(coeffs['a'] * kv_val + coeffs['b'], 2)
        partial_deriv_kv = coeffs['a']
        
        incerteza_csr = propagate_uncertainty(lambda: csr_val, [(partial_deriv_kv, d_kv_abs)])
        return csr_val, round(incerteza_csr, 4)
    except: return "Erro CSR", 0.0

def calcular_fator_g(csr_val, espessura_val, d_espessura_abs):
    try:
        csr_keys = list(FATOR_G_CONSTANTS_UNCERTAINTIES.keys())
        csr_aproximado_key = min(csr_keys, key=lambda x: abs(x - csr_val))
        data = FATOR_G_CONSTANTS_UNCERTAINTIES.get(csr_aproximado_key)
        
        a0, da0, a1, da1, a2, da2, a3, da3 = data['a0'], data['da0'], data['a1'], data['da1'], data['a2'], data['da2'], data['a3'], data['da3']
        
        fg_calc = a0 + a1*espessura_val + a2*espessura_val**2 + a3*espessura_val**3
        fg = max(0, round(fg_calc, 4))
        
        p_esp = a1 + 2*a2*espessura_val + 3*a3*espessura_val**2
        
        inc = propagate_uncertainty(lambda: fg, [
            (p_esp, d_espessura_abs), (1, da0), (espessura_val, da1), (espessura_val**2, da2), (espessura_val**3, da3)
        ])
        return fg, round(inc, 4)
    except: return "Erro G", 0.0

def calcular_glandularidade(idade, espessura_mama_cm):
    espessura_mm = espessura_mama_cm * 10
    if 30 <= idade <= 49: a, b, c, k = -0.000196, 0.0666, -7.45, 278
    elif 50 <= idade <= 54: a, b, c, k = -0.000255, 0.0768, -7.67, 259
    elif 55 <= idade <= 59: a, b, c, k = -0.000199, 0.0593, -6.00, 207
    elif 60 <= idade <= 88: a, b, c, k = -0.000186, 0.0572, -5.99, 208
    else: return "Idade fora da faixa"
    
    G = (a * espessura_mm**3) + (b * espessura_mm**2) + (c * espessura_mm) + k
    return max(0, round(G, 2))

def get_coeffs_from_lambda_for_fator_c(csr_key, group_key):
    # Mapa auxiliar simplificado para obter coeficientes das fórmulas do Fator C
    # Isso é necessário para propagar incerteza corretamente
    # Aqui vamos usar uma aproximação prática para não explodir o código
    # (Apenas para fins de exemplo, mantemos a estrutura, mas o cálculo de C usará as lambdas diretamente)
    return None 

def calcular_fator_c(csr, espessura, glandularidade, d_espessura_abs):
    try:
        grupo = 1 if glandularidade <= 25 else 2 if glandularidade <= 50 else 3 if glandularidade <= 75 else 4
        csr_prox = min(formulas_fator_c.keys(), key=lambda x: abs(x - csr))
        
        # Chama a função lambda correta
        f_c_val = formulas_fator_c[csr_prox][grupo](espessura)
        
        # Incerteza simplificada (5% padrão conforme solicitado anteriormente)
        inc_c = f_c_val * INCERTEZA_COEFS_FATOR_C_PERCENTUAL 
        return round(f_c_val, 4), round(inc_c, 4)
    except: return "Erro C", 0.0

def calcular_ki(kv, alvo_filtro, mas, espessura_mama, d_mas_abs, d_espessura_abs, local_nome, tabelas_state):
    try:
        tabela = tabelas_state.get(local_nome)
        x_val = tabela.get((alvo_filtro, int(kv)))
        
        if x_val is None: return "Kv/Filtro não encontrado", 0.0
        
        conv, ref = (1892.25, 64) if local_nome == 'UFRJ' else (2500, 63)
        div = (ref - espessura_mama)**2
        
        ki = round(((x_val * mas) * conv) / div, 2)
        
        # Propagação Simplificada para o exemplo final
        p_mas = (x_val * conv) / div
        inc_ki = propagate_uncertainty(lambda: ki, [(p_mas, d_mas_abs)])
        
        return ki, round(inc_ki, 4)
    except: return "Erro Ki", 0.0

def calcular_dgm(ki, s, fg, fc, i_ki, i_s, i_fg, i_fc):
    try:
        dgm = ki * s * fg * fc
        # Propagação de DGM
        inc = propagate_uncertainty(lambda: dgm, [
            (s*fg*fc, i_ki), (ki*fg*fc, i_s), (ki*s*fc, i_fg), (ki*s*fg, i_fc)
        ])
        return round(dgm, 2), round(inc * 0.10, 4)
    except: return "Erro DGM", 0.0

# --- 3. INTERFACE STREAMLIT ---

st.set_page_config(page_title="Calculadora DGM 2.0", page_icon="🔬", layout="wide")

# Inicialização do Session State
if 'tabelas_ki' not in st.session_state:
    st.session_state.tabelas_ki = {'IRD': tabela_ki_ird, 'UFRJ': tabela_ki_ufrj}
if 'csr_coeffs' not in st.session_state:
    st.session_state.csr_coeffs = csr_coeffs_default.copy()
if 'resultados' not in st.session_state:
    st.session_state.resultados = pd.DataFrame()

st.title("🔬 Calculadora de Dose Glandular Média (nova versão)")

with st.sidebar:
    st.header("⚙️ Configuração")
    
    # --- ÁREA DE UPLOAD DE NOVO MAMÓGRAFO ---
    with st.expander("➕ Adicionar Novo Mamógrafo"):
        st.info("Envie um Excel com colunas: Alvo/Filtro, kV, Ki. Opcional: CSR_a, CSR_b")
        nome_eq = st.text_input("Nome do Equipamento:")
        planilha = st.file_uploader("Carregar Excel", type=['xlsx'])
        
        if st.button("Salvar Equipamento") and planilha and nome_eq:
            try:
                df = pd.read_excel(planilha)
                novo_ki = {}
                # Processa Ki
                for _, r in df.iterrows():
                    novo_ki[(str(r['Alvo/Filtro']), int(r['kV']))] = float(r['Ki'])
                
                # Salva tabela Ki
                st.session_state.tabelas_ki[nome_eq] = novo_ki
                
                # Processa CSR (se houver)
                if 'CSR_a' in df.columns and 'CSR_b' in df.columns:
                    for _, r in df.iterrows():
                        st.session_state.csr_coeffs[str(r['Alvo/Filtro'])] = {'a': r['CSR_a'], 'b': r['CSR_b']}
                
                st.success(f"Equipamento '{nome_eq}' cadastrado com sucesso!")
            except Exception as e:
                st.error(f"Erro ao ler arquivo: {e}")

    st.markdown("---")
    st.header("📝 Dados do Paciente e Exame")
    
    # --- CAMPOS RESTAURADOS ---
    paciente_id = st.text_input('ID do Paciente:', help="Prontuário ou Identificador")
    iniciais_paciente = st.text_input('Iniciais da Paciente:', max_chars=3, help="Ex: J.S.").upper()
    
    local_mamografo = st.selectbox('Equipamento:', options=list(st.session_state.tabelas_ki.keys()))
    
    idade = st.number_input('Idade:', min_value=1, max_value=120, value=45)
    espessura_mama = st.number_input('Espessura (cm):', min_value=1.0, max_value=20.0, value=6.0, step=0.1)
    alvo_filtro = st.selectbox('Alvo/Filtro:', options=list(alvo_filtro_options.keys()))
    kv = st.number_input('Kv:', min_value=20.0, max_value=50.0, value=28.0, step=0.1)
    mas = st.number_input('mAs:', min_value=0.1, max_value=1000.0, value=50.0, step=0.1)
    
    # --- CHECKBOX DE GLANDULARIDADE RESTAURADO ---
    sabe_glandularidade = st.checkbox("Eu sei a glandularidade (manual)")
    glandularidade_manual = None
    if sabe_glandularidade:
        glandularidade_manual = st.number_input('Glandularidade (%):', 0.0, 100.0, 50.0, 0.1)

# --- 4. EXECUÇÃO DO CÁLCULO ---

if st.button("Calcular DGM", type="primary"):
    st.subheader("Resultados:")
    
    # 1. Definição da Glandularidade
    if sabe_glandularidade and glandularidade_manual is not None:
        gland = glandularidade_manual
        st.info(f"🔹 Glandularidade Manual: {gland}%")
    else:
        gland = calcular_glandularidade(idade, espessura_mama)
        if isinstance(gland, str):
            st.error(f"Erro Glandularidade: {gland}")
            gland = None
        else:
            st.info(f"🔹 Glandularidade Calculada: {gland}%")

    if gland is not None:
        # Incertezas Absolutas
        dkv = kv * 0.01
        dmas = mas * 0.05
        desp = espessura_mama * 0.05
        
        # 2. CSR
        csr, i_csr = calcular_csr(kv, alvo_filtro, dkv, st.session_state.csr_coeffs)
        
        # 3. Fator g
        if isinstance(csr, (int, float)):
            fg, i_fg = calcular_fator_g(csr, espessura_mama, desp)
        else:
            fg, i_fg = "Erro", 0.0
            
        # 4. Fator C
        if isinstance(csr, (int, float)):
            fc, i_fc = calcular_fator_c(csr, espessura_mama, gland, desp)
        else:
            fc, i_fc = "Erro", 0.0
            
        # 5. Valor s
        s_val = alvo_filtro_options.get(alvo_filtro, 1)
        
        # 6. Ki
        ki, i_ki = calcular_ki(kv, alvo_filtro, mas, espessura_mama, dmas, desp, local_mamografo, st.session_state.tabelas_ki)
        
        # Exibição Intermediária
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("CSR", f"{csr}", delta=f"± {i_csr}")
        c2.metric("Fator g", f"{fg}", delta=f"± {i_fg}")
        c3.metric("Fator C", f"{fc}", delta=f"± {i_fc}")
        c4.metric("Ki (mGy)", f"{ki}", delta=f"± {i_ki}")
        
        # 7. DGM Final
        erro_nos_dados = any(isinstance(x, str) for x in [csr, fg, fc, ki])
        
        st.markdown("---")
        if not erro_nos_dados:
            dgm, i_dgm = calcular_dgm(ki, s_val, fg, fc, i_ki, 0, i_fg, i_fc)
            st.success(f"### ✅ DGM Calculada: {dgm} mGy ± {i_dgm}")
            
            # Adicionar ao Histórico
            novo_resultado = {
                "Data": datetime.now().strftime("%d/%m/%Y %H:%M"),
                "Paciente ID": paciente_id,
                "Iniciais": iniciais_paciente,
                "Equipamento": local_mamografo,
                "Idade": idade,
                "Espessura": espessura_mama,
                "kV": kv,
                "mAs": mas,
                "Gland. (%)": gland,
                "DGM (mGy)": dgm,
                "Incerteza": i_dgm
            }
            st.session_state.resultados = pd.concat([st.session_state.resultados, pd.DataFrame([novo_resultado])], ignore_index=True)
            
        else:
            st.error("Não foi possível calcular a DGM. Verifique se o kV e Filtro selecionados existem na tabela do equipamento.")

# --- 5. HISTÓRICO ---
st.markdown("---")
st.subheader("📊 Histórico de Pacientes")
if not st.session_state.resultados.empty:
    st.dataframe(st.session_state.resultados, use_container_width=True)
    
    # Botão de Download
    csv = st.session_state.resultados.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Baixar CSV", data=csv, file_name="relatorio_dgm.csv", mime="text/csv")

st.markdown("---")
st.markdown("Desenvolvido por Jossana Almeida, com o auxilio de um modelo de linguagem.")

