"""
Que hace este programa?

Ejecuta la comparación entre:
  - Algoritmo GENÉTICO ORIGINAL
  - Algoritmo GENÉTICO MEJORADO (elitismo + torneo + mutación inteligente)
  
Características principales:
  - La versión original usa penalización por gap = -4.
  - La versión mejorada usa penalización por gap = -2.
  - Ambas series se normalizan con el mismo offset para graficar en positivo.
  - Se valida la integridad de las secuencias (sin gaps) al final.
Requisitos: módulo 'blosum' (BLOSUM(62)), matplotlib.
"""

import random
import copy
import time
import matplotlib.pyplot as plt
import blosum

GENERACIONES = 100
N_INDIVIDUOS_ORIG = 10    # población inicial en el original (con mutación inicial)
N_INDIVIDUOS_MEJ = 20     # población inicial en el mejorado (copias de originales)
N_GAPS_INICIAL = 1        # gaps insertados en la inicial del original (tal como tu código)
PROB_MUT_INDIVIDUO = 0.8  # probabilidad de mutar en la versión original (igual que tu código)
TORNEO_K = 3              # tamaño del torneo para selección en versión mejorada
PORC_SELECCION = 0.5      # porcentaje que se mantiene tras eliminar_peores (truncation)
ELITISMO = True           # activar elitismo en la versión mejorada

blosum62 = blosum.BLOSUM(62)
NFE = 0

def get_sequences():
    seq1 = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"
    seq2 = "MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"
    seq3 = "MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"
    return [list(seq1), list(seq2), list(seq3)]


def crear_individuo():
    return [row[:] for row in get_sequences()]

def crear_poblacion_inicial_original(n=10):
    individuo_base = crear_individuo()
    return [[row[:] for row in individuo_base] for _ in range(n)]

def igualar_longitud_secuencias(individuo, gap='-'):
    max_len = max(len(fila) for fila in individuo)
    return [fila + [gap] * (max_len - len(fila)) for fila in individuo]

def mutar_poblacion_v2(poblacion, num_gaps=1):
    poblacion_mutada = []
    for individuo in poblacion:
        nuevo_individuo = []
        for fila in individuo:
            fila_mutada = fila[:]
            posiciones = set()
            for _ in range(num_gaps):
                pos = random.randint(0, len(fila_mutada))
                while pos in posiciones:
                    pos = random.randint(0, len(fila_mutada))
                posiciones.add(pos)
                fila_mutada.insert(pos, '-')
            nuevo_individuo.append(fila_mutada)
        poblacion_mutada.append(nuevo_individuo)
    return poblacion_mutada

def eliminar_peores(poblacion, scores, porcentaje=0.5):
    idx_ordenados = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    n_seleccionados = max(1, int(len(poblacion) * porcentaje))
    ind_seleccionados = [poblacion[i] for i in idx_ordenados[:n_seleccionados]]
    scores_seleccionados = [scores[i] for i in idx_ordenados[:n_seleccionados]]
    return ind_seleccionados, scores_seleccionados

def mutar_individuo(individuo, n_gaps, p):
    nuevo_individuo = []
    for secuencia in individuo:
        sec = secuencia[:]
        if random.random() < p:
            posiciones = set()
            for _ in range(n_gaps):
                pos = random.randint(0, len(sec))
                while pos in posiciones:
                    pos = random.randint(0, len(sec))
                posiciones.add(pos)
                sec.insert(pos, '-')
        nuevo_individuo.append(sec)
    return nuevo_individuo

def cruzar_individuos_doble_punto(ind1, ind2):
    hijo1 = []
    hijo2 = []
    for seq1, seq2 in zip(ind1, ind2):
        aa_indices = [i for i, a in enumerate(seq1) if a != '-']
        if len(aa_indices) < 6:
            hijo1.append(seq1[:])
            hijo2.append(seq2[:])
            continue
        intentos = 0
        while True:
            p1, p2 = sorted(random.sample(aa_indices, 2))
            if p2 - p1 >= 5 or intentos > 10:
                break
            intentos += 1

        def cruza(seqA, seqB):
            aaA = [a for a in seqA if a != '-']
            aaB = [a for a in seqB if a != '-']
            nueva = aaA[:p1] + aaB[p1:p2] + aaA[p2:]
            resultado = []
            idx = 0
            for a in seqA:
                if a == '-':
                    resultado.append('-')
                else:
                    resultado.append(nueva[idx])
                    idx += 1
            return resultado

        nueva_seq1 = cruza(seq1, seq2)
        nueva_seq2 = cruza(seq2, seq1)
        hijo1.append(nueva_seq1)
        hijo2.append(nueva_seq2)

    hijo1 = mutar_individuo(hijo1, 1, 0.8)
    hijo2 = mutar_individuo(hijo2, 1, 0.8)
    return hijo1, hijo2

def cruzar_poblacion_doble_punto(poblacion):
    nueva_poblacion = []
    n = len(poblacion)
    indices = list(range(n))
    random.shuffle(indices)
    parejas = [(indices[i], indices[i+1]) for i in range(0, n-1, 2)]
    if n % 2 == 1:
        parejas.append((indices[-1], indices[0]))
    for idx1, idx2 in parejas:
        padre1 = poblacion[idx1]
        padre2 = poblacion[idx2]
        hijo1, hijo2 = cruzar_individuos_doble_punto(padre1, padre2)
        nueva_poblacion.append(copy.deepcopy(padre1))
        nueva_poblacion.append(copy.deepcopy(padre2))
        nueva_poblacion.append(hijo1)
        nueva_poblacion.append(hijo2)
    return nueva_poblacion[:2*n]

def obtener_best(scores, poblacion):
    idx_mejor = scores.index(max(scores))
    fitness_best = scores[idx_mejor]
    best = copy.deepcopy(poblacion[idx_mejor])
    return best, fitness_best

def validar_poblacion_sin_gaps(poblacion, originales):
    for individuo in poblacion:
        for seq, seq_orig in zip(individuo, originales):
            seq_sin_gaps = [a for a in seq if a != '-']
            seq_orig_sin_gaps = [a for a in seq_orig if a != '-']
            if seq_sin_gaps != seq_orig_sin_gaps:
                return False
    return True


def evaluar_individuo_blosum(individuo, gap_penalty=-4):
    """
    Evaluación general parametrizable:
      - gap_penalty: penalización por gap (usamos -4 en el original, -2 en el mejorado)
    Devuelve el score (puede ser negativo).
    """
    global NFE
    NFE += 1
    score = 0
    n_seqs = len(individuo)
    seq_len = len(individuo[0])
    for col in range(seq_len):
        for i in range(n_seqs):
            for j in range(i+1, n_seqs):
                a = individuo[i][col]
                b = individuo[j][col]
                if a == '-' or b == '-':
                    score += gap_penalty
                else:
                    score += blosum62[a][b]
    return score

def ejecutar_ag_original(n_individuos=10, num_gaps=1, prob_mut=0.8, n_generaciones=100):
    """
    Ejecuta exactamente la versión original (con mutación inicial).
    gap_penalty en la evaluación se fija a -4 aquí, que era tu original.
    """
    # 1) crear población inicial y mutar como lo hacías
    poblacion = crear_poblacion_inicial_original(n_individuos)
    poblacion = mutar_poblacion_v2(poblacion, num_gaps=num_gaps)
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]

    # 2) evaluar y seleccionar
    scores = [evaluar_individuo_blosum(ind, gap_penalty=-4) for ind in poblacion]
    poblacion, scores = eliminar_peores(poblacion, scores, porcentaje=PORC_SELECCION)

    fitness_por_generacion = []
    veryBest = None
    fitnessVeryBest = None

    for gen in range(n_generaciones):
        poblacion = cruzar_poblacion_doble_punto(poblacion)
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
        scores = [evaluar_individuo_blosum(ind, gap_penalty=-4) for ind in poblacion]
        poblacion, scores = eliminar_peores(poblacion, scores, porcentaje=PORC_SELECCION)
        best, fitness_best = obtener_best(scores, poblacion)
        if veryBest is None or fitness_best > fitnessVeryBest:
            veryBest = best
            fitnessVeryBest = fitness_best
        fitness_por_generacion.append(fitnessVeryBest)

    return fitness_por_generacion, veryBest


def mutacion_inteligente_evaluada(individuo, gap_penalty, prob=0.4):
    # intenta eliminar o insertar gaps; acepta cambios solo si no empeoran el score
    base_copia = [fila[:] for fila in individuo]
    base_copia = igualar_longitud_secuencias(base_copia)
    base_score = evaluar_individuo_blosum(base_copia, gap_penalty=gap_penalty)

    nuevo = [fila[:] for fila in individuo]

    for s_idx, seq in enumerate(individuo):
        # intento de eliminación de gap
        gap_positions = [i for i, a in enumerate(seq) if a == '-']
        if gap_positions and random.random() < prob:
            pos = random.choice(gap_positions)
            candidato = [fila[:] for fila in nuevo]
            candidato[s_idx].pop(pos)
            candidato = igualar_longitud_secuencias(candidato)
            sc = evaluar_individuo_blosum(candidato, gap_penalty=gap_penalty)
            if sc >= base_score:
                nuevo = candidato
                base_score = sc

        # intento de insertar gap (prob pequeña)
        if random.random() < (prob * 0.5):
            pos = random.randint(0, len(seq))
            candidato = [fila[:] for fila in nuevo]
            candidato[s_idx].insert(pos, '-')
            candidato = igualar_longitud_secuencias(candidato)
            sc = evaluar_individuo_blosum(candidato, gap_penalty=gap_penalty)
            if sc >= base_score:
                nuevo = candidato
                base_score = sc

    # asegurar que lo que retornamos tiene longitudes iguales
    nuevo = igualar_longitud_secuencias(nuevo)
    return nuevo


def ejecutar_ag_mejorado(n_individuos=20, n_gaps_mut=1, prob_mut=0.8, n_generaciones=100,
                         gap_penalty=-2, torneo_k=3, porcentaje_sel=0.5, elitismo=True, verbose=False):
    """
    - Poblacion inicial: copias de los 3 originales (no mutar al inicio)
    - Selección: torneo para elegir padres (presión selectiva)
    - Cruza: se reutiliza cruzar_individuos_doble_punto
    - Mutación: mutacion_inteligente_evaluada (elimina gaps dañinos y solo acepta cambios que mejoren o igualen)
    - gap_penalty: penalización por gap (aquí la opción A: -2)
    """
    # 1) población inicial (copias de los originales, sin cambios)
    base = crear_individuo()
    poblacion = [copy.deepcopy(base) for _ in range(n_individuos)]
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]

    # eval inicial y truncation selection
    scores = [evaluar_individuo_blosum(ind, gap_penalty=gap_penalty) for ind in poblacion]
    poblacion, scores = eliminar_peores(poblacion, scores, porcentaje=porcentaje_sel)

    fitness_por_generacion = []
    veryBest = None
    fitnessVeryBest = None
    start_time = time.time()

    for gen in range(n_generaciones):
        # 1) crear pool de padres por torneo
        padres = []
        for _ in range(len(poblacion)):
            idxs = random.sample(range(len(poblacion)), min(torneo_k, len(poblacion)))
            ganador_idx = max(idxs, key=lambda idx: scores[idx])
            padres.append(copy.deepcopy(poblacion[ganador_idx]))

        # 2) cruza por parejas (aleatorizar orden)
        random.shuffle(padres)
        poblacion_cruce = []
        for i in range(0, len(padres)-1, 2):
            p1 = padres[i]
            p2 = padres[i+1]
            h1, h2 = cruzar_individuos_doble_punto(p1, p2)
            # aplicar mutación inteligente evaluada con prob adaptativa (reduce con gen)
            prob_gen = prob_mut * (1 - gen / max(1, n_generaciones))
            h1 = mutacion_inteligente_evaluada(h1, gap_penalty=gap_penalty, prob=prob_gen)
            h2 = mutacion_inteligente_evaluada(h2, gap_penalty=gap_penalty, prob=prob_gen)
            poblacion_cruce.extend([copy.deepcopy(p1), copy.deepcopy(p2), h1, h2])

        # si impar, emparejar último con uno aleatorio
        if len(padres) % 2 == 1:
            p_last = padres[-1]
            p_rand = random.choice(padres[:-1])
            h1, h2 = cruzar_individuos_doble_punto(p_last, p_rand)
            prob_gen = prob_mut * (1 - gen / max(1, n_generaciones))
            h1 = mutacion_inteligente_evaluada(h1, gap_penalty=gap_penalty, prob=prob_gen)
            h2 = mutacion_inteligente_evaluada(h2, gap_penalty=gap_penalty, prob=prob_gen)
            poblacion_cruce.extend([copy.deepcopy(p_last), copy.deepcopy(p_rand), h1, h2])

        # 3) igualar longitudes y evaluar
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion_cruce]
        scores = [evaluar_individuo_blosum(ind, gap_penalty=gap_penalty) for ind in poblacion]

        # 4) selección truncada
        poblacion, scores = eliminar_peores(poblacion, scores, porcentaje=porcentaje_sel)

        # 5) elitismo: si historial es mejor que actuales, reinyectar
        best, fitness_best = obtener_best(scores, poblacion)
        if veryBest is None or fitness_best > fitnessVeryBest:
            veryBest = best
            fitnessVeryBest = fitness_best
        else:
            if elitismo and fitnessVeryBest > max(scores):
                # reinsertar elite en la última posición
                poblacion[-1] = copy.deepcopy(veryBest)
                scores[-1] = fitnessVeryBest

        fitness_por_generacion.append(fitnessVeryBest)
        if verbose and (gen % 10 == 0 or gen == n_generaciones-1):
            elapsed = time.time() - start_time
            print(f"[Mejorado] Gen {gen+1}/{n_generaciones} best={fitnessVeryBest} time={elapsed:.2f}s")

    return fitness_por_generacion, veryBest

if __name__ == "__main__":
    print("Ejecutando algoritmo ORIGINAL (penal gap = -4)...")
    fitness_orig_raw, best_orig = ejecutar_ag_original(n_individuos=N_INDIVIDUOS_ORIG,
                                                       num_gaps=N_GAPS_INICIAL,
                                                       prob_mut=PROB_MUT_INDIVIDUO,
                                                       n_generaciones=GENERACIONES)

    print("\nEjecutando algoritmo MEJORADO (penal gap = -2, elitismo, torneo, mutación inteligente)...")
    fitness_mej_raw, best_mej = ejecutar_ag_mejorado(n_individuos=N_INDIVIDUOS_MEJ,
                                                     n_gaps_mut=N_GAPS_INICIAL,
                                                     prob_mut=PROB_MUT_INDIVIDUO,
                                                     n_generaciones=GENERACIONES,
                                                     gap_penalty=-2,
                                                     torneo_k=TORNEO_K,
                                                     porcentaje_sel=PORC_SELECCION,
                                                     elitismo=ELITISMO,
                                                     verbose=True)

    all_values = []
    all_values.extend(fitness_orig_raw)
    all_values.extend(fitness_mej_raw)
    min_val = min(all_values) if all_values else 0
    offset = -min_val + 1 if min_val <= 0 else 0

    fitness_orig_pos = [f + offset for f in fitness_orig_raw]
    fitness_mej_pos = [f + offset for f in fitness_mej_raw]

  
    plt.figure(figsize=(10,6))
    plt.plot(range(len(fitness_orig_pos)), fitness_orig_pos, label='Original (penal -4)')
    plt.plot(range(len(fitness_mej_pos)), fitness_mej_pos, label='Mejorado (penal -2)')
    plt.xlabel('Generaciones')
    plt.ylabel('Fitness (normalizado a positivo, mayor = mejor)')
    plt.title('Comparación: Algoritmo Original vs Algoritmo Mejorado')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


    originales = crear_individuo()
    def sin_gaps_a_strs(ind):
        return ["".join(seq).replace('-', '') for seq in ind]

    best_orig_nogap = sin_gaps_a_strs(best_orig)
    best_mej_nogap = sin_gaps_a_strs(best_mej)
    originales_strs = ["".join(seq) for seq in originales]

    valid_orig = all(o == b for o, b in zip(originales_strs, best_orig_nogap))
    valid_mej = all(o == b for o, b in zip(originales_strs, best_mej_nogap))

    print("\nValidación de integridad:")
    print(f" - Mejor ORIGINAL conserva secuencias base (sin gaps eliminados)? {valid_orig}")
    print(f" - Mejor MEJORADO conserva secuencias base (sin gaps eliminados)? {valid_mej}")

    print("\nResumen final:")
    print(f" - Último fitness raw (original): {fitness_orig_raw[-1] if fitness_orig_raw else None}")
    print(f" - Último fitness raw (mejorado): {fitness_mej_raw[-1] if fitness_mej_raw else None}")
    print(f" - Offset aplicado para graficar en positivo: {offset}")


