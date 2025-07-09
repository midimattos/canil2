import streamlit as st
import sqlite3
from collections import Counter
from itertools import product
import json

# --- Configuração SQLite ---
conn = sqlite3.connect("caes.db", check_same_thread=False)
cursor = conn.cursor()
cursor.execute("""
CREATE TABLE IF NOT EXISTS caes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    nome TEXT UNIQUE,
    pelagem TEXT,
    trufagem TEXT,
    coloracao TEXT,
    marcacao TEXT,
    diluicao TEXT,
    intensidade TEXT,
    olhos TEXT,
    genotipo TEXT
)
""")
conn.commit()

# --- Funções de genética ---

def inferir_genotipo(dados):
    genotipo = {}
    coloracao = dados['coloracao']
    marcacao = dados['marcacao']
    diluicao = dados['diluicao']
    olhos = dados['olhos']

    genotipo['Locus_E'] = 'E/E' if coloracao != 'Laranja/Creme' else 'e/e'
    genotipo['Locus_K'] = 'K/K' if coloracao in ['Preto Sólido', 'Chocolate'] else 'k/k'

    if coloracao == 'Preto Sólido':
        genotipo['Locus_A'] = 'a/a'
        genotipo['Locus_B'] = 'B/B'
    elif coloracao == 'Chocolate':
        genotipo['Locus_A'] = 'a/a'
        genotipo['Locus_B'] = 'b/b'
    elif coloracao == 'Laranja/Creme':
        genotipo['Locus_A'] = 'ay/ay'
        genotipo['Locus_B'] = 'B/B'
    elif coloracao == 'Agouti':
        genotipo['Locus_A'] = 'aw/aw'
        genotipo['Locus_B'] = 'B/B'
    else:
        genotipo['Locus_A'] = 'a/a'
        genotipo['Locus_B'] = 'B/B'

    genotipo['Locus_D'] = 'd/d' if diluicao == 'Sim' else 'D/D'

    if marcacao == 'Particolor':
        genotipo['Locus_S'] = 'sp/sp'
    elif marcacao == 'Tan Points':
        genotipo['Locus_S'] = 'at/at'
    elif marcacao == 'Agouti':
        genotipo['Locus_S'] = 'aw/aw'
    else:
        genotipo['Locus_S'] = 'S/S'

    genotipo['Locus_M'] = 'm/m'
    genotipo['Locus_I'] = 'i/i' if olhos == 'Azuis' else 'I/I'
    return genotipo

def genotipo_para_caracteristicas(genotipo):
    caracteristicas = {}
    caracteristicas['Pelagem'] = 'Longa'

    B = genotipo['Locus_B']
    D = genotipo['Locus_D']

    if 'b/b' == B and D == 'd/d':
        caracteristicas['Trufagem'] = 'Marrom Claro (Diluição)'
    elif 'b/b' == B:
        caracteristicas['Trufagem'] = 'Marrom'
    else:
        caracteristicas['Trufagem'] = 'Preta'

    A = genotipo['Locus_A']
    if A == 'ay/ay':
        cor = 'Laranja/Creme'
    elif A == 'aw/aw':
        cor = 'Agouti'
    elif B == 'b/b' and D == 'd/d':
        cor = 'Chocolate Diluído'
    elif B == 'b/b':
        cor = 'Chocolate'
    else:
        cor = 'Preto Sólido'
    caracteristicas['Coloração'] = cor

    S = genotipo['Locus_S']
    if S == 'sp/sp':
        marcacao = 'Particolor'
    elif S == 'at/at':
        marcacao = 'Tan Points'
    elif S == 'aw/aw':
        marcacao = 'Agouti'
    else:
        marcacao = 'Nenhuma'
    caracteristicas['Marcação'] = marcacao

    caracteristicas['Diluição'] = 'Sim' if D == 'd/d' else 'Não'

    I = genotipo['Locus_I']
    if I == 'i/i':
        caracteristicas['Olhos'] = 'Azuis'
    else:
        caracteristicas['Olhos'] = 'Escuros'

    caracteristicas['Intensidade'] = 'Baixa' if caracteristicas['Olhos'] == 'Azuis' else 'Alta'
    return caracteristicas

def punnett_combinations(ale1, ale2):
    return [f"{min(a1,a2)}/{max(a1,a2)}" for a1 in ale1 for a2 in ale2]

def calcular_probabilidades(pai_gen, mae_gen):
    loci = pai_gen.keys()
    combinacoes_locus = {}
    for locus in loci:
        alelos_pai = pai_gen[locus].split('/')
        alelos_mae = mae_gen[locus].split('/')
        combos = punnett_combinations(alelos_pai, alelos_mae)
        freq = Counter(combos)
        total = sum(freq.values())
        probabilidade = {k: v/total for k,v in freq.items()}
        combinacoes_locus[locus] = probabilidade

    loci_list = list(loci)
    todas_combinacoes = list(product(*[list(combinacoes_locus[l].items()) for l in loci_list]))

    resultados = []
    for comb in todas_combinacoes:
        genotipo = {}
        prob = 1
        for i, (alelos, p) in enumerate(comb):
            genotipo[loci_list[i]] = alelos
            prob *= p

        # 🚫 FILTROS DE IMPOSSIBILIDADE GENÉTICA
        B_ = genotipo['Locus_B']
        D_ = genotipo['Locus_D']
        I_ = genotipo['Locus_I']

        # Chocolate inviável se não existir b/b
        if B_ == 'b/b' and ('b' not in pai_gen['Locus_B'] and 'b' not in mae_gen['Locus_B']):
            continue
        # Diluição inviável se ambos forem D/D
        if D_ == 'd/d' and (pai_gen['Locus_D'] == 'D/D' and mae_gen['Locus_D'] == 'D/D'):
            continue
        # Olhos azuis inviável se ambos forem I/I
        if I_ == 'i/i' and (pai_gen['Locus_I'] == 'I/I' and mae_gen['Locus_I'] == 'I/I'):
            continue

        resultados.append((genotipo, prob))

    resultados.sort(key=lambda x: x[1], reverse=True)
    return resultados

def agrupar_por_caracteristicas(resultados):
    grupos = {}
    for gen, prob in resultados:
        carac = tuple(sorted(genotipo_para_caracteristicas(gen).items()))
        grupos[carac] = grupos.get(carac, 0) + prob
    grupos_ordenados = sorted(grupos.items(), key=lambda x: x[1], reverse=True)
    return grupos_ordenados

# --- CSS para esconder sidebar no mobile e ajustar área principal ---
st.markdown(
    """
    <style>
    /* Esconde sidebar para telas menores que 600px */
    @media (max-width: 600px) {
        .css-18e3th9 {
            display: none;
        }
        .css-1d391kg {
            width: 100% !important;
            padding-left: 1rem !important;
            padding-right: 1rem !important;
        }
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# --- Menu horizontal no topo para todos ---
menu = st.radio("Menu", ["Cadastrar Cão", "Catálogo", "Simular Ninhada"], horizontal=True)

if menu == "Cadastrar Cão":
    st.header("📋 Cadastro de Cão")
    nome = st.text_input("Nome do cão")
    pelagem = st.selectbox("Pelagem", ["Longa"])
    trufagem = st.selectbox("Trufagem", ["Preta", "Marrom"])
    coloracao = st.selectbox("Coloração", ["Preto Sólido", "Chocolate", "Laranja/Creme", "Agouti"])
    marcacao = st.selectbox("Marcação", ["Nenhuma", "Tan Points", "Agouti", "Particolor"])
    diluicao = st.selectbox("Diluição", ["Não", "Sim"])
    intensidade = st.selectbox("Intensidade", ["Alta", "Baixa"])
    olhos = st.selectbox("Olhos", ["Escuros", "Azuis"])

    if st.button("Cadastrar"):
        dados = {
            'pelagem': pelagem,
            'trufagem': trufagem,
            'coloracao': coloracao,
            'marcacao': marcacao,
            'diluicao': diluicao,
            'intensidade': intensidade,
            'olhos': olhos
        }
        genotipo = inferir_genotipo(dados)
        genotipo_json = json.dumps(genotipo)
        try:
            cursor.execute(
                "INSERT INTO caes (nome, pelagem, trufagem, coloracao, marcacao, diluicao, intensidade, olhos, genotipo) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (nome, pelagem, trufagem, coloracao, marcacao, diluicao, intensidade, olhos, genotipo_json))
            conn.commit()
            st.success(f"Cão {nome} cadastrado com sucesso!")
        except sqlite3.IntegrityError:
            st.error("Já existe um cão cadastrado com esse nome.")

elif menu == "Catálogo":
    st.header("🐾 Catálogo de Cães")
    cursor.execute("SELECT nome, pelagem, trufagem, coloracao, marcacao, diluicao, intensidade, olhos FROM caes ORDER BY nome")
    caes = cursor.fetchall()
    if not caes:
        st.info("Nenhum cão cadastrado ainda.")
    else:
        for c in caes:
            with st.expander(c[0]):
                st.write(f"Pelagem: {c[1]}")
                st.write(f"Trufagem: {c[2]}")
                st.write(f"Coloração: {c[3]}")
                st.write(f"Marcação: {c[4]}")
                st.write(f"Diluição: {c[5]}")
                st.write(f"Intensidade: {c[6]}")
                st.write(f"Olhos: {c[7]}")

elif menu == "Simular Ninhada":
    st.header("🔮 Simulador de Ninhada")
    cursor.execute("SELECT nome, genotipo FROM caes ORDER BY nome")
    caes = cursor.fetchall()

    if len(caes) < 2:
        st.warning("Cadastre pelo menos dois cães para simular uma ninhada.")
    else:
        # Criar dicionário seguro para busca de genótipos
        dicionario_caes = {nome: genotipo for nome, genotipo in caes}

        nomes_caes = list(dicionario_caes.keys())
        pai_nome = st.selectbox("Selecione o Pai", nomes_caes)
        mae_nome = st.selectbox("Selecione a Mãe", nomes_caes)

        if pai_nome == mae_nome:
            st.warning("Pai e mãe não podem ser o mesmo cão.")
        elif st.button("Simular"):
            pai_gen = json.loads(dicionario_caes[pai_nome])
            mae_gen = json.loads(dicionario_caes[mae_nome])

            resultados = calcular_probabilidades(pai_gen, mae_gen)
            grupos = agrupar_por_caracteristicas(resultados)

            st.subheader(f"🐶 Possíveis filhotes entre {pai_nome} e {mae_nome}:")

            if not grupos:
                st.error("❌ Nenhum filhote viável pode nascer dessa combinação, de acordo com a genética.")
            else:
                # Mostrar no máximo 3 filhotes mais prováveis
                num_filhotes = min(3, len(grupos))
                for i, (carac, prob) in enumerate(grupos[:num_filhotes], start=1):
                    st.write(f"### 🍼 Filhote {i} - Chance: {prob*100:.2f}%")
                    for k, v in carac:
                        st.write(f"**{k}:** {v}")
                    st.write("---")
