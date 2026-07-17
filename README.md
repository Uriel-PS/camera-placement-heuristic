# Surveillance System: Heuristic Camera Placement at Critical Points

## 1. The Problem

This project implements a variation of the classic **Art Gallery Problem**, which is **NP-hard**: given a polygonal environment (here represented as a matrix/grid), the goal is to determine the minimum number of "guards" (cameras) — and where to place them — to cover the entire visible area of the environment.

There is no known polynomial-time algorithm that solves this problem exactly for large instances (the search space is combinatorial: for a matrix with N free cells and 4 possible camera directions, the number of possible subsets is on the order of `2^(4N)`). For this reason, a **greedy constructive heuristic** combined with **local search** is applied here, which produces good solutions in polynomial time, without guaranteeing optimality.

## 2. Real-World Variables Incorporated

| # | Real-World Variable | How It Was Modeled |
|---|---|---|
| A | **Line of sight / obstacles** | Each grid cell can be `0` (free) or `1` (wall/pillar). The algorithm uses **raycasting** (Bresenham's algorithm) to trace a straight line between the camera and each candidate cell; if the line crosses a wall before reaching the target, the cell is in **shadow** (blind spot), even if it's within the camera's nominal range. |
| B | **Limited field of view (FOV)** | Real cameras don't instantly rotate 360°. Each camera has a **fixed direction** (North, South, East, West) and a **configurable opening angle** (e.g. 90°). A cell is only seen if it's within the FOV cone AND within maximum range (radius, in cells) AND has an unobstructed line of sight (variable A). |

These two variables make the problem far closer to reality than a simple "Wi-Fi signal" passing through walls or a 360° camera with unlimited range.

## 3. Modeling

- **Grid**: `rows x columns` matrix, `0` = free space, `1` = wall.
- **Camera candidate**: pair `(free position, direction ∈ {N,S,E,W})`.
- **Coverage of a camera**: set of cells within range R, within the FOV cone, and with unobstructed line of sight (no walls in the path).
- **Objective**: given a maximum number of cameras `K`, maximize the total number of covered cells (union of all cameras' coverage).

This is essentially a **Set Cover** problem, which is NP-hard, with the peculiarity that the "sets" (visible areas) change shape depending on the environment's geometry (obstacles) and the chosen direction, unlike classic Set Cover, where sets are fixed.

## 4. Heuristic Algorithm

### Phase 1 — Greedy Construction (Greedy Set Cover)
At each iteration:
1. For **every** free cell and **each** of the 4 directions, calculate how many **currently uncovered** cells that camera would cover (marginal gain).
2. Choose the (position, direction) combination with the **highest marginal gain**.
3. Mark those cells as covered and repeat until reaching the maximum number of cameras, 100% coverage, or no further gain is possible.

This is the classic heuristic for Set Cover, guaranteeing an approximation factor of `ln(n)` relative to the optimal solution in the general case.

### Phase 2 — Local Search (refinement)
After the greedy phase, the algorithm attempts to escape simple local minima: for each already-placed camera, it is **temporarily removed** and the algorithm searches for the **best alternative position/direction** considering the coverage of the other fixed cameras. If a different position increases total coverage, the camera is **relocated**. The process repeats until convergence (no camera changes position) or a maximum number of iterations is reached.

### Complexity
- Exact solution (brute force): exponential, `O(2^(4N))`.
- Greedy heuristic: approximately `O(K · N · R²)`, where `K` = number of cameras, `N` = grid cells, `R` = range — **polynomial**.
- Local search: adds a small multiplicative factor (few iterations).

## 5. Code Structure (`vigilancia.c`)

- `linha_de_visao_livre`: Raycasting (Bresenham) — variable (A).
- `dentro_do_fov`: checks whether a cell is within the FOV cone — variable (B).
- `calcula_visibilidade`: combines range + FOV + raycasting to determine everything a hypothetical camera would see.
- `heuristica_gulosa`: Phase 1.
- `busca_local`: Phase 2.
- `imprime_mapa` / `imprime_relatorio`: ASCII visualization and final report.
- `carrega_mapa_arquivo` / `carrega_mapa_exemplo`: data input.

## 6. Building and Running

```bash
gcc -O2 -Wall -o vigilancia vigilancia.c -lm

# Using the built-in example map (warehouse with pillars):
./vigilancia

# Using a custom map, with custom parameters:
./vigilancia mapa_banco.txt <max_cameras> <fov_degrees> <range>
./vigilancia mapa_banco.txt 6 90 6
```

### Map File Format
Plain text, `0` = free, `1` = wall, space-separated values, one file line = one matrix row (see the included `mapa_banco.txt`, which represents a small bank branch with rooms and a hallway).

## 7. Example Output (summary)

With the example map (15x20 warehouse, 148 free cells), 12 cameras, 90° FOV, range 8: the algorithm achieves **74.3%** coverage, leaving 38 blind spots (cells behind pillars/corners that would require more cameras or wider FOV). On the bank branch map (7x12, 35 free cells), 6 cameras with range 6 achieve **94.3%** coverage.

The final report lists the exact position and direction of each camera, and the ASCII map shows:
- `#` wall
- `.` covered cell
- `?` blind spot (not covered)
- `^ > v <` camera position and direction (North/East/South/West)

This allows direct visualization of remaining blind spots and evaluation of whether they're acceptable or justify purchasing additional cameras — a real cost-benefit decision in a physical security project.

## Authors

Uriel Pacheco de Souza, [Kelwin Efrain Bagnhuk da Silva](https://www.linkedin.com/in/kelwin-efrain-bagnhuk-da-silva-4a9331407/), and Emanuel

## Documentation

The full academic report (in Portuguese) is available in this repository as a PDF.
