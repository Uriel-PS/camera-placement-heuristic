# Sistema de Vigilância: Alocação Heurística de Câmeras em Pontos Críticos

## 1. O Problema

Este trabalho implementa uma variação do clássico **Problema da Galeria de Arte**
(*Art Gallery Problem*), que é **NP-difícil**: dado um ambiente poligonal (aqui
representado como uma matriz/grid), deseja-se determinar o menor número de
"guardas" (câmeras) — e onde posicioná-los — de forma a cobrir toda a área
visível do ambiente.

Não existe algoritmo conhecido de tempo polinomial que resolva esse problema de
forma exata para instâncias grandes (o espaço de busca é combinatório: para uma
matriz com N células livres e 4 direções possíveis de câmera, o número de
subconjuntos possíveis é da ordem de `2^(4N)`). Por isso, aqui é aplicada uma
**heurística construtiva gulosa** combinada com **busca local**, que produz boas
soluções em tempo polinomial, mas sem garantia de otimalidade.

## 2. Variáveis do mundo real incorporadas

| # | Variável real | Como foi modelada |
|---|----------------|--------------------|
| A | **Linha de visão / obstáculos** | Cada célula da matriz pode ser `0` (livre) ou `1` (parede/pilastra). O algoritmo usa **Raycasting** (algoritmo de Bresenham) para traçar uma linha reta entre a câmera e cada célula candidata; se a linha cruzar uma parede antes de chegar ao destino, a célula fica em **sombra** (ponto cego), mesmo estando dentro do alcance nominal da câmera. |
| B | **Campo de visão limitado (FOV)** | Câmeras reais não giram 360° instantaneamente. Cada câmera tem uma **direção fixa** (Norte, Sul, Leste, Oeste) e um **ângulo de abertura configurável** (ex: 90°). Uma célula só é vista se estiver dentro do cone de abertura E dentro do alcance máximo (raio, em número de células) E com linha de visão livre (variável A). |

Essas duas variáveis tornam o problema muito mais próximo da realidade do que um
simples "Wi-Fi" que atravessa paredes ou uma câmera 360° sem limite de alcance —
exatamente como pedido no enunciado (cf. exemplo do Caixeiro Viajante com grafo
não completo e assimetria de rotas).

## 3. Modelagem

- **Grid**: matriz `linhas x colunas`, `0` = espaço livre, `1` = parede.
- **Candidato de câmera**: par `(posição livre, direção ∈ {N,S,L,O})`.
- **Cobertura de uma câmera**: conjunto de células dentro do alcance R, dentro
  do cone de FOV, e com linha de visão desobstruída (sem paredes no caminho).
- **Objetivo**: dado um número máximo de câmeras `K`, maximizar o número total
  de células cobertas (união das coberturas de todas as câmeras).

Este é essencialmente um problema de **Cobertura de Conjuntos (Set Cover)**, que
é NP-difícil, com a peculiaridade de que os "conjuntos" (áreas visíveis) mudam
de formato dependendo da geometria do ambiente (obstáculos) e da direção
escolhida — diferente do Set Cover clássico, onde os conjuntos são fixos.

## 4. Algoritmo Heurístico

### Fase 1 — Construção Gulosa (Greedy Set Cover)
A cada iteração:
1. Para **cada** célula livre e **cada** uma das 4 direções, calcula-se quantas
   células **ainda não cobertas** aquela câmera cobriria (ganho marginal).
2. Escolhe-se a combinação (posição, direção) de **maior ganho marginal**.
3. Marca essas células como cobertas e repete até atingir o número máximo de
   câmeras ou 100% de cobertura (ou nenhum ganho adicional ser possível).

Esta é a heurística clássica para Set Cover, que garante um fator de
aproximação de `ln(n)` em relação à solução ótima no caso geral.

### Fase 2 — Busca Local (refinamento)
Após a fase gulosa, o algoritmo tenta escapar de mínimos locais simples:
para cada câmera já posicionada, ela é **removida temporariamente** e o
algoritmo busca a **melhor posição/direção alternativa** considerando a
cobertura das demais câmeras fixas. Se uma posição diferente aumentar a
cobertura total, a câmera é **realocada**. O processo se repete até convergir
(nenhuma câmera muda de posição) ou atingir um número máximo de iterações.

### Complexidade
- Solução exata (força bruta): exponencial, `O(2^(4N))`.
- Heurística gulosa: aproximadamente `O(K · N · R²)`, onde `K` = número de
  câmeras, `N` = células do grid, `R` = alcance — **polinomial**.
- Busca local: adiciona um fator multiplicativo pequeno (poucas iterações).

## 5. Estrutura do código (`vigilancia.c`)

- `linha_de_visao_livre`: Raycasting (Bresenham) — variável (A).
- `dentro_do_fov`: verifica ângulo dentro do cone de visão — variável (B).
- `calcula_visibilidade`: combina alcance + FOV + raycasting para determinar
  tudo que uma câmera hipotética enxergaria.
- `heuristica_gulosa`: Fase 1.
- `busca_local`: Fase 2.
- `imprime_mapa` / `imprime_relatorio`: visualização ASCII e relatório final.
- `carrega_mapa_arquivo` / `carrega_mapa_exemplo`: entrada de dados.

## 6. Compilação e execução

```bash
gcc -O2 -Wall -o vigilancia vigilancia.c -lm

# Usando o mapa de exemplo embutido (galpão com pilastras):
./vigilancia

# Usando um mapa próprio, com parâmetros customizados:
./vigilancia mapa_banco.txt <max_cameras> <fov_graus> <alcance>
./vigilancia mapa_banco.txt 6 90 6
```

### Formato do arquivo de mapa
Texto simples, `0` = livre, `1` = parede, valores separados por espaço, uma
linha do arquivo = uma linha da matriz (veja `mapa_banco.txt` incluso, que
representa uma pequena agência bancária com salas e um corredor).

## 7. Exemplo de saída (resumo)

Com o mapa de exemplo (galpão 15x20, 148 células livres), 12 câmeras, FOV 90°,
alcance 8: o algoritmo atinge **74,3%** de cobertura, deixando 38 pontos cegos
(células atrás de pilastras/cantos que exigiriam mais câmeras ou câmeras com
FOV maior). Já no mapa da agência bancária (7x12, 35 células livres), 6
câmeras com alcance 6 atingem **94,3%** de cobertura.

O relatório final lista a posição e direção exata de cada câmera, e o mapa
ASCII mostra:
- `#` parede
- `.` célula coberta
- `?` ponto cego (não coberto)
- `^ > v <` posição e direção da câmera (Norte/Leste/Sul/Oeste)

Isso permite visualizar diretamente onde estão os pontos cegos remanescentes
e avaliar se são aceitáveis ou se justificam a compra de câmeras adicionais —
uma decisão real de custo-benefício em um projeto de segurança física.
