/* ============================================================================
 * SISTEMA DE VIGILANCIA: Alocacao Heuristica de Cameras em Pontos Criticos
 * ============================================================================
 *
 * PROBLEMA:
 *   Variacao do Problema da Galeria de Arte (Art Gallery Problem), que e
 *   NP-Dificil. Dado um galpao/planta representado por uma matriz com
 *   paredes (obstaculos), queremos posicionar um numero limitado de cameras
 *   (ou o menor numero possivel) de forma a maximizar a area coberta,
 *   respeitando:
 *
 *   (A) LINHA DE VISAO REAL: a camera nao enxerga atraves de paredes. Usamos
 *       Raycasting (algoritmo de Bresenham) para tracar raios da camera ate
 *       cada celula candidata; se o raio colide com uma parede antes de
 *       chegar la, a celula fica em "sombra" (ponto cego).
 *
 *   (B) CAMPO DE VISAO LIMITADO (FOV): cameras reais nao giram 360 graus.
 *       Cada camera tem uma direcao fixa (N/S/L/O) e um angulo de abertura
 *       (ex: 90 graus) e um alcance maximo (raio). So enxerga o que estiver
 *       dentro desse cone/leque, alem de respeitar obstaculos.
 *
 * ABORDAGEM (heuristica, pois o problema exato e NP-Difícil):
 *   1. Geracao de candidatos: toda celula livre e um possivel local de
 *      camera, em qualquer uma das 4 direcoes.
 *   2. Heuristica GULOSA (Greedy Set-Cover): a cada passo, escolhe a
 *      combinacao (posicao, direcao) que cobre o MAIOR NUMERO DE CELULAS
 *      AINDA NAO COBERTAS. Repete ate atingir cobertura desejada ou o
 *      numero maximo de cameras.
 *   3. Busca Local (melhoria 2-opt): apos a fase gulosa, tenta reposicionar
 *      cada camera individualmente (posicao e direcao) para uma configuracao
 *      que aumente a cobertura total, escapando de minimos locais simples
 *      do guloso.
 *
 * COMPLEXIDADE:
 *   Testar todas as combinacoes (subconjuntos de posicoes x direcoes) seria
 *   exponencial -> O(2^N). A heuristica gulosa roda em tempo polinomial
 *   (aproximadamente O(K * N * R^2), K=numero de cameras, N=celulas,
 *   R=alcance), entregando uma solucao boa, mas nao necessariamente otima.
 *
 * COMPILACAO:
 *   gcc -O2 -Wall -o vigilancia vigilancia.c -lm
 *
 * EXECUCAO:
 *   ./vigilancia [arquivo_mapa.txt] [max_cameras] [fov_graus] [alcance]
 *   ./vigilancia                      -> usa mapa de exemplo embutido
 *
 * FORMATO DO MAPA (arquivo texto):
 *   0 = espaco livre    1 = parede/obstaculo
 *   Cada linha do arquivo = uma linha da matriz, numeros separados por espaco
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_DIM 60
#define MAX_CAMERAS_LIMIT 200
#define PI 3.14159265358979323846

typedef enum { NORTE = 0, LESTE = 1, SUL = 2, OESTE = 3, NUM_DIRECOES = 4 } Direcao;

const char *NOME_DIRECAO[NUM_DIRECOES] = {"Norte", "Leste", "Sul", "Oeste"};
/* angulo central de cada direcao, em graus (0 = Leste, sentido anti-horario) */
const double ANGULO_CENTRAL[NUM_DIRECOES] = {90.0, 0.0, 270.0, 180.0};

typedef struct {
    int linha, coluna;
    Direcao direcao;
    int celulas_cobertas; /* quantas celulas NOVAS esta camera cobriu quando foi colocada */
} Camera;

typedef struct {
    int grid[MAX_DIM][MAX_DIM];      /* 0 = livre, 1 = parede */
    int linhas, colunas;
    int coberto[MAX_DIM][MAX_DIM];   /* 0/1 - celula esta vigiada por alguma camera */
    int fov_graus;                   /* angulo de abertura da camera, ex 90 */
    int alcance;                     /* raio maximo de visao em celulas */
    Camera cameras[MAX_CAMERAS_LIMIT];
    int num_cameras;
    int total_livres;                /* total de celulas livres no mapa */
} Mapa;

/* ---------------------------------------------------------------------
 * Utilitarios de angulo
 * --------------------------------------------------------------------- */

/* Normaliza um angulo para o intervalo [0, 360) */
double normaliza_angulo(double a) {
    while (a < 0) a += 360.0;
    while (a >= 360.0) a -= 360.0;
    return a;
}

/* Verifica se o angulo 'alvo' esta dentro do cone de FOV centrado em 'centro' */
int dentro_do_fov(double centro, double alvo, double fov) {
    double diff = normaliza_angulo(alvo - centro);
    if (diff > 180.0) diff = 360.0 - diff;
    return diff <= (fov / 2.0 + 1e-9);
}

/* ---------------------------------------------------------------------
 * RAYCASTING (Bresenham) - traca uma linha entre dois pontos e verifica
 * se ha colisao com parede antes de alcancar o destino.
 * Retorna 1 se a linha de visao esta LIVRE (sem obstaculos no caminho),
 * 0 se algo bloqueia.
 * --------------------------------------------------------------------- */
int linha_de_visao_livre(Mapa *m, int r0, int c0, int r1, int c1) {
    int dr = abs(r1 - r0), dc = abs(c1 - c0);
    int sr = (r0 < r1) ? 1 : -1;
    int sc = (c0 < c1) ? 1 : -1;
    int err = dr - dc;
    int r = r0, c = c0;

    while (r != r1 || c != c1) {
        int e2 = 2 * err;
        if (e2 > -dc) { err -= dc; r += sr; }
        if (e2 <  dr) { err += dr; c += sc; }

        if (r == r1 && c == c1) break;

        /* celula intermediaria bloqueada -> sem visao */
        if (m->grid[r][c] == 1) return 0;
    }
    return 1;
}

/* ---------------------------------------------------------------------
 * Calcula o conjunto de celulas visiveis por UMA camera hipotetica em
 * (linha, coluna) apontando para 'dir', respeitando FOV, alcance e paredes.
 * Preenche o vetor 'visivel' (mesma dimensao do grid) com 1 nas celulas
 * visiveis. Retorna a quantidade de celulas visiveis.
 * --------------------------------------------------------------------- */
int calcula_visibilidade(Mapa *m, int linha, int coluna, Direcao dir,
                          int visivel[MAX_DIM][MAX_DIM]) {
    memset(visivel, 0, sizeof(int) * MAX_DIM * MAX_DIM);
    if (m->grid[linha][coluna] == 1) return 0; /* camera nao fica dentro de parede */

    double centro = ANGULO_CENTRAL[dir];
    int count = 0;

    int rmin = linha - m->alcance, rmax = linha + m->alcance;
    int cmin = coluna - m->alcance, cmax = coluna + m->alcance;
    if (rmin < 0) rmin = 0;
    if (cmin < 0) cmin = 0;
    if (rmax >= m->linhas) rmax = m->linhas - 1;
    if (cmax >= m->colunas) cmax = m->colunas - 1;

    for (int r = rmin; r <= rmax; r++) {
        for (int c = cmin; c <= cmax; c++) {
            if (m->grid[r][c] == 1) continue;          /* celula-alvo e parede */
            int dr = r - linha, dc = c - coluna;
            double dist = sqrt((double)(dr * dr + dc * dc));
            if (dist > m->alcance) continue;             /* fora do alcance */

            if (dr == 0 && dc == 0) {
                /* a propria posicao da camera esta sempre visivel */
                visivel[r][c] = 1;
                count++;
                continue;
            }

            /* angulo do alvo em relacao a camera (eixo Y invertido: linha
               cresce para "baixo" no grid, entao usamos -dr para Norte ficar
               para cima visualmente) */
            double angulo_alvo = normaliza_angulo(atan2(-dr, dc) * 180.0 / PI);

            if (!dentro_do_fov(centro, angulo_alvo, m->fov_graus)) continue;

            if (linha_de_visao_livre(m, linha, coluna, r, c)) {
                visivel[r][c] = 1;
                count++;
            }
        }
    }
    return count;
}

/* ---------------------------------------------------------------------
 * Carrega mapa de um arquivo texto. Retorna 1 em sucesso, 0 em erro.
 * --------------------------------------------------------------------- */
int carrega_mapa_arquivo(Mapa *m, const char *caminho) {
    FILE *f = fopen(caminho, "r");
    if (!f) return 0;

    char linha_buf[1024];
    int r = 0;
    while (fgets(linha_buf, sizeof(linha_buf), f) && r < MAX_DIM) {
        int c = 0;
        char *tok = strtok(linha_buf, " \t\r\n");
        while (tok && c < MAX_DIM) {
            m->grid[r][c] = atoi(tok);
            c++;
            tok = strtok(NULL, " \t\r\n");
        }
        if (c == 0) continue; /* linha vazia, ignora */
        m->colunas = c;
        r++;
    }
    m->linhas = r;
    fclose(f);
    return 1;
}

/* Mapa de exemplo embutido, representando um galpao/banco com pilastras
   (paredes internas) e uma sala separada - simula um cenario real onde
   corredores tem cantos e obstaculos que criam pontos cegos. */
void carrega_mapa_exemplo(Mapa *m) {
    int exemplo[15][20] = {
        {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1},
        {1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,0,1},
        {1,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1},
        {1,0,1,0,1,1,1,0,1,0,0,1,0,1,1,0,1,0,1,1},
        {1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1},
        {1,1,1,0,1,0,1,1,1,1,1,1,1,0,1,0,1,1,0,1},
        {1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1},
        {1,0,1,1,1,0,1,0,1,1,0,1,1,0,1,1,0,1,0,1},
        {1,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1},
        {1,1,1,0,1,1,1,0,1,0,1,1,1,1,1,0,1,1,0,1},
        {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
        {1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1},
        {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
    };
    m->linhas = 15;
    m->colunas = 20;
    for (int r = 0; r < m->linhas; r++)
        for (int c = 0; c < m->colunas; c++)
            m->grid[r][c] = exemplo[r][c];
}

void inicializa_mapa(Mapa *m, int fov_graus, int alcance) {
    m->fov_graus = fov_graus;
    m->alcance = alcance;
    m->num_cameras = 0;
    memset(m->coberto, 0, sizeof(m->coberto));
    m->total_livres = 0;
    for (int r = 0; r < m->linhas; r++)
        for (int c = 0; c < m->colunas; c++)
            if (m->grid[r][c] == 0) m->total_livres++;
}

int total_coberto(Mapa *m) {
    int total = 0;
    for (int r = 0; r < m->linhas; r++)
        for (int c = 0; c < m->colunas; c++)
            if (m->grid[r][c] == 0 && m->coberto[r][c]) total++;
    return total;
}

/* Aplica a cobertura de uma camera ao mapa->coberto (marca as celulas) */
void aplica_cobertura(Mapa *m, int visivel[MAX_DIM][MAX_DIM]) {
    for (int r = 0; r < m->linhas; r++)
        for (int c = 0; c < m->colunas; c++)
            if (visivel[r][c]) m->coberto[r][c] = 1;
}

/* ---------------------------------------------------------------------
 * HEURISTICA GULOSA (Greedy Set Cover)
 * A cada iteracao, testa TODAS as posicoes livres x 4 direcoes e escolhe
 * a que cobre o maior numero de celulas AINDA NAO cobertas (ganho marginal
 * maximo). Repete ate atingir max_cameras ou 100% de cobertura.
 * --------------------------------------------------------------------- */
void heuristica_gulosa(Mapa *m, int max_cameras, double cobertura_alvo) {
    int visivel_tmp[MAX_DIM][MAX_DIM];

    while (m->num_cameras < max_cameras) {
        int melhor_ganho = 0;
        int melhor_r = -1, melhor_c = -1;
        Direcao melhor_dir = NORTE;

        for (int r = 0; r < m->linhas; r++) {
            for (int c = 0; c < m->colunas; c++) {
                if (m->grid[r][c] == 1) continue; /* nao coloca camera em parede */

                for (int d = 0; d < NUM_DIRECOES; d++) {
                    calcula_visibilidade(m, r, c, (Direcao)d, visivel_tmp);

                    /* calcula ganho = celulas visiveis que AINDA NAO estavam cobertas */
                    int ganho = 0;
                    for (int rr = 0; rr < m->linhas; rr++)
                        for (int cc = 0; cc < m->colunas; cc++)
                            if (visivel_tmp[rr][cc] && !m->coberto[rr][cc]) ganho++;

                    if (ganho > melhor_ganho) {
                        melhor_ganho = ganho;
                        melhor_r = r; melhor_c = c; melhor_dir = (Direcao)d;
                    }
                }
            }
        }

        if (melhor_ganho <= 0) {
            printf("[Guloso] Nenhuma posicao adicional melhora a cobertura. Parando.\n");
            break;
        }

        /* efetiva a melhor camera encontrada */
        calcula_visibilidade(m, melhor_r, melhor_c, melhor_dir, visivel_tmp);
        aplica_cobertura(m, visivel_tmp);

        m->cameras[m->num_cameras].linha = melhor_r;
        m->cameras[m->num_cameras].coluna = melhor_c;
        m->cameras[m->num_cameras].direcao = melhor_dir;
        m->cameras[m->num_cameras].celulas_cobertas = melhor_ganho;
        m->num_cameras++;

        double pct = 100.0 * total_coberto(m) / m->total_livres;
        printf("  Camera %2d colocada em (%2d,%2d) direcao %-6s | +%3d celulas | cobertura total: %.1f%%\n",
               m->num_cameras, melhor_r, melhor_c, NOME_DIRECAO[melhor_dir], melhor_ganho, pct);

        if (pct >= cobertura_alvo) {
            printf("[Guloso] Cobertura alvo de %.1f%% atingida.\n", cobertura_alvo);
            break;
        }
    }
}

/* ---------------------------------------------------------------------
 * Recalcula toda a matriz 'coberto' a partir do zero, considerando o
 * conjunto atual de cameras. Necessario apos mover/remover uma camera
 * durante a busca local.
 * --------------------------------------------------------------------- */
void recalcula_cobertura_total(Mapa *m) {
    int visivel_tmp[MAX_DIM][MAX_DIM];
    memset(m->coberto, 0, sizeof(m->coberto));
    for (int i = 0; i < m->num_cameras; i++) {
        calcula_visibilidade(m, m->cameras[i].linha, m->cameras[i].coluna,
                              m->cameras[i].direcao, visivel_tmp);
        aplica_cobertura(m, visivel_tmp);
    }
}

/* ---------------------------------------------------------------------
 * BUSCA LOCAL (melhoria pos-guloso, estilo 2-opt):
 * Para cada camera ja posicionada, tenta REMOVE-LA temporariamente e
 * procurar a MELHOR posicao/direcao alternativa (considerando a cobertura
 * das demais cameras fixas). Se a alternativa gerar cobertura total maior
 * que a atual, a camera e realocada. Repete por algumas rodadas ou ate
 * nao haver mais melhoria (convergencia local).
 * --------------------------------------------------------------------- */
void busca_local(Mapa *m, int max_iteracoes) {
    int visivel_tmp[MAX_DIM][MAX_DIM];
    printf("\n[Busca Local] Iniciando refinamento das posicoes...\n");

    for (int iter = 0; iter < max_iteracoes; iter++) {
        int houve_melhoria = 0;

        for (int i = 0; i < m->num_cameras; i++) {
            Camera original = m->cameras[i];

            /* remove temporariamente a camera i (substitui por sentinela) */
            Camera guardada = m->cameras[i];
            /* recalcula cobertura das OUTRAS cameras (sem a i) */
            memset(m->coberto, 0, sizeof(m->coberto));
            for (int j = 0; j < m->num_cameras; j++) {
                if (j == i) continue;
                calcula_visibilidade(m, m->cameras[j].linha, m->cameras[j].coluna,
                                      m->cameras[j].direcao, visivel_tmp);
                aplica_cobertura(m, visivel_tmp);
            }
            int cobertura_sem_i = total_coberto(m);

            /* procura a melhor posicao/direcao para reinserir a camera i */
            int melhor_ganho = -1, melhor_r = original.linha, melhor_c = original.coluna;
            Direcao melhor_dir = original.direcao;

            for (int r = 0; r < m->linhas; r++) {
                for (int c = 0; c < m->colunas; c++) {
                    if (m->grid[r][c] == 1) continue;
                    for (int d = 0; d < NUM_DIRECOES; d++) {
                        calcula_visibilidade(m, r, c, (Direcao)d, visivel_tmp);
                        int ganho = 0;
                        for (int rr = 0; rr < m->linhas; rr++)
                            for (int cc = 0; cc < m->colunas; cc++)
                                if (visivel_tmp[rr][cc] && !m->coberto[rr][cc]) ganho++;
                        if (ganho > melhor_ganho) {
                            melhor_ganho = ganho;
                            melhor_r = r; melhor_c = c; melhor_dir = (Direcao)d;
                        }
                    }
                }
            }

            /* aplica a melhor camera encontrada (pode ser a mesma de antes) */
            m->cameras[i].linha = melhor_r;
            m->cameras[i].coluna = melhor_c;
            m->cameras[i].direcao = melhor_dir;
            calcula_visibilidade(m, melhor_r, melhor_c, melhor_dir, visivel_tmp);
            aplica_cobertura(m, visivel_tmp);

            int cobertura_nova_total = cobertura_sem_i + melhor_ganho;

            if (melhor_r != guardada.linha || melhor_c != guardada.coluna ||
                melhor_dir != guardada.direcao) {
                houve_melhoria = 1;
                printf("  [iter %d] Camera %d realocada: (%d,%d,%s) -> (%d,%d,%s) | cobertura total: %d/%d\n",
                       iter + 1, i + 1,
                       guardada.linha, guardada.coluna, NOME_DIRECAO[guardada.direcao],
                       melhor_r, melhor_c, NOME_DIRECAO[melhor_dir],
                       cobertura_nova_total, m->total_livres);
            }
        }

        recalcula_cobertura_total(m);

        if (!houve_melhoria) {
            printf("[Busca Local] Convergiu (sem melhorias) na iteracao %d.\n", iter + 1);
            break;
        }
    }
}

/* ---------------------------------------------------------------------
 * VISUALIZACAO EM ASCII
 *   #  = parede
 *   .  = espaco coberto por alguma camera
 *   ?  = espaco livre NAO coberto (ponto cego)
 *   0-9/A.. = posicao da camera (numero de identificacao) + seta de direcao
 * --------------------------------------------------------------------- */
void imprime_mapa(Mapa *m) {
    char tela[MAX_DIM][MAX_DIM];
    for (int r = 0; r < m->linhas; r++) {
        for (int c = 0; c < m->colunas; c++) {
            if (m->grid[r][c] == 1) tela[r][c] = '#';
            else if (m->coberto[r][c]) tela[r][c] = '.';
            else tela[r][c] = '?';
        }
    }
    /* marca as cameras por cima, com seta indicando direcao */
    const char setas[NUM_DIRECOES] = {'^', '>', 'v', '<'}; /* N, L, S, O */
    for (int i = 0; i < m->num_cameras; i++) {
        int r = m->cameras[i].linha, c = m->cameras[i].coluna;
        tela[r][c] = setas[m->cameras[i].direcao];
    }

    printf("\n     ");
    for (int c = 0; c < m->colunas; c++) printf("%d", c % 10);
    printf("\n");
    for (int r = 0; r < m->linhas; r++) {
        printf("%3d  ", r);
        for (int c = 0; c < m->colunas; c++) putchar(tela[r][c]);
        printf("\n");
    }
    printf("\nLegenda: # parede | . coberto | ? ponto cego (nao coberto) | ^ > v < camera (direcao)\n");
}

void imprime_relatorio(Mapa *m) {
    int coberto = total_coberto(m);
    printf("\n================= RELATORIO FINAL =================\n");
    printf("Dimensao do mapa: %d x %d\n", m->linhas, m->colunas);
    printf("Celulas livres (area util): %d\n", m->total_livres);
    printf("Cameras utilizadas: %d\n", m->num_cameras);
    printf("FOV por camera: %d graus | Alcance: %d celulas\n", m->fov_graus, m->alcance);
    printf("Celulas cobertas: %d\n", coberto);
    printf("Pontos cegos remanescentes: %d\n", m->total_livres - coberto);
    printf("Cobertura total: %.2f%%\n", 100.0 * coberto / m->total_livres);
    printf("=====================================================\n");
    printf("\nPosicoes finais das cameras:\n");
    for (int i = 0; i < m->num_cameras; i++) {
        printf("  Camera %2d: linha=%2d coluna=%2d direcao=%-6s\n",
               i + 1, m->cameras[i].linha, m->cameras[i].coluna,
               NOME_DIRECAO[m->cameras[i].direcao]);
    }
}

/* ---------------------------------------------------------------------
 * MAIN
 * --------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    Mapa mapa;
    int max_cameras = 12;
    int fov_graus = 90;
    int alcance = 8;
    double cobertura_alvo = 100.0; /* tenta cobrir 100%, para quando nao ha mais ganho */

    if (argc >= 2) {
        if (!carrega_mapa_arquivo(&mapa, argv[1])) {
            fprintf(stderr, "Erro ao abrir arquivo de mapa '%s'. Usando mapa de exemplo.\n", argv[1]);
            carrega_mapa_exemplo(&mapa);
        }
    } else {
        carrega_mapa_exemplo(&mapa);
    }
    if (argc >= 3) max_cameras = atoi(argv[2]);
    if (argc >= 4) fov_graus = atoi(argv[3]);
    if (argc >= 5) alcance = atoi(argv[4]);

    if (max_cameras > MAX_CAMERAS_LIMIT) max_cameras = MAX_CAMERAS_LIMIT;

    inicializa_mapa(&mapa, fov_graus, alcance);

    printf("============================================================\n");
    printf(" SISTEMA DE VIGILANCIA - Alocacao Heuristica de Cameras\n");
    printf("============================================================\n");
    printf("Mapa: %d x %d | Celulas livres: %d\n", mapa.linhas, mapa.colunas, mapa.total_livres);
    printf("Parametros: max_cameras=%d | FOV=%d graus | alcance=%d celulas\n\n",
           max_cameras, fov_graus, alcance);

    printf("Mapa original (paredes # / espaco livre ?):\n");
    imprime_mapa(&mapa);

    printf("\n--- FASE 1: Heuristica Gulosa (Greedy Set Cover) ---\n");
    heuristica_gulosa(&mapa, max_cameras, cobertura_alvo);

    printf("\nMapa apos fase gulosa:\n");
    imprime_mapa(&mapa);

    printf("\n--- FASE 2: Busca Local (refinamento) ---\n");
    busca_local(&mapa, 5);

    printf("\nMapa apos busca local (resultado final):\n");
    imprime_mapa(&mapa);

    imprime_relatorio(&mapa);

    return 0;
}
