# ARQUITETURA-MODULAR-DE-FUSAO-DE-SENSORES-PARA-NAVEGACAO-AUTONOMA-SEGURANCA-EFICIENCIA-E-ETICA
# Simulação de fusão de sensores com EKF e monitor de consistência

Código de apoio ao artigo **"Arquitetura modular de fusão de sensores para navegação autônoma: segurança, eficiência e ética"** (CIITS'26, Universidade CEUMA).

Simula um veículo planar com IMU, odometria e GNSS e compara quatro métodos de estimação de posição em três cenários, usando repetições Monte Carlo com sementes fixas. Todos os dados são **sintéticos**.

## Conteúdo

| Arquivo | Função |
|---|---|
| `simulacao_ekf.py` | Simulação principal: gera a trajetória, os sensores, roda os quatro métodos nos três cenários e produz `resultados.json` e as figuras |
| `sensibilidade.py` | Análise de sensibilidade: varia o ruído real do GNSS (2,5; 4,0; 6,0 m) mantendo em 2,5 m o valor assumido pelo filtro. Gera `sensibilidade.json` |
| `resultados.json` | Resultados agregados da simulação principal (média, desvio-padrão e mediana de 100 execuções) |
| `sensibilidade.json` | Resultados da análise de sensibilidade |

## Requisitos

- Python 3.9 ou superior
- `numpy`
- `matplotlib`

```bash
pip install numpy matplotlib
```

## Como executar

Execute os dois scripts na mesma pasta, pois `sensibilidade.py` importa funções de `simulacao_ekf.py`.

```bash
python simulacao_ekf.py     # gera resultados.json e as figuras (.png)
python sensibilidade.py     # gera sensibilidade.json e imprime um resumo
```

A execução leva alguns minutos, dependendo da máquina (100 execuções × 3 cenários × 2 versões do EKF, mais 300 execuções na sensibilidade).

### Saídas

- `resultados.json`: métricas por cenário e método.
- `fig_trajetoria.png`: trajetória verdadeira, fixes do GNSS, estima por odometria e EKF (cenário de queda do GNSS, 1ª execução).
- `fig_erro_tempo.png`: erro de posição ao longo do tempo e raio de incerteza de 95 % (queda do GNSS, 1ª execução).
- `fig_rmse.png`: RMSE por método e cenário, em escala logarítmica.
- `sensibilidade.json`: RMSE, cobertura e taxa de rejeição para cada ruído real do GNSS.

## O que é simulado

### Trajetória

Veículo planar por 300 s, velocidade de 8 ± 1,5 m/s e taxa de giro senoidal composta. O estimador roda a 20 Hz (passo de 0,05 s).

### Sensores

| Sensor | Taxa | Modelo |
|---|---|---|
| GNSS | 1 Hz | Posição + ruído gaussiano (σ = 2,5 m por eixo) |
| Odometria | 10 Hz | Velocidade × (1 + erro de escala de 0,5 %) + ruído (σ = 0,15 m/s) |
| Giroscópio | 20 Hz | Taxa de giro + viés constante por execução (σ = 0,005 rad/s) + ruído (σ = 0,01 rad/s) |
| Acelerômetro | 20 Hz | Aceleração + ruído (σ = 0,10 m/s²) |

### Cenários

| Cenário | Descrição |
|---|---|
| `nominal` | Sem falhas |
| `queda_gnss` | GNSS indisponível de 100 s a 140 s |
| `outliers_gnss` | 8 % dos fixos recebem deslocamento aleatório de 25 a 40 m (o primeiro fixo é sempre válido) |

### Métodos comparados

1. **GNSS isolado**: mantém o último fixo recebido.
2. **Estima por odometria**: integra giroscópio e velocidade das rodas, sem correção absoluta, inicializada no primeiro fixo.
3. **EKF sem validação**: EKF com estado `[px, py, ψ, v, b_g]` (posição, orientação, velocidade, viés do giroscópio), atualizações de GNSS e odometria com covariância na forma de Joseph.
4. **EKF com validação (NIS)**: como o anterior, mas rejeita fixos do GNSS cuja inovação normalizada ao quadrado excede 9,21 (qui-quadrado, 2 gl, 99 %). Após 5 rejeições consecutivas, reinicializa a posição com a medição corrente (recuperação de travamento).

### Métricas

| Chave no JSON | Significado |
|---|---|
| `rmse`, `emax` | Erro quadrático médio e erro máximo de posição (m) |
| `path_err_pct` | Erro relativo do comprimento da trajetória estimada (%). Indicador indireto de correções de controle; **não é consumo de energia** |
| `cobertura95` | % do tempo em que o erro real não excede o raio de incerteza de 95 % (ideal: ≈ 95 %) |
| `pct_degradado` | % do tempo com r95 > 5 m (modo degradado) |
| `emax_sem_alerta` | Erro máximo enquanto o sistema está em modo normal |
| `latencia_alerta_s` | Tempo entre o início da queda do GNSS e o alerta (só no cenário `queda_gnss`) |
| `rej_pct` | % de fixos do GNSS rejeitados pela validação |
| `t_cycle_us` | Tempo de processamento por ciclo do EKF (µs) |

Cada métrica é agregada em `mean`, `std` e `median` sobre as 100 execuções.

## Parâmetros principais

Ficam no início de `simulacao_ekf.py`:

| Constante | Valor | Descrição |
|---|---|---|
| `DT` | 0,05 s | Passo do estimador |
| `T_TOTAL` | 300 s | Duração |
| `SIG_GPS` | 2,5 m | Desvio-padrão do GNSS (assumido pelo filtro) |
| `GATE` | 9,21 | Limiar de validação (NIS) |
| `R95_LIM` | 5,0 m | Limiar do modo degradado |
| `N_RUNS` | 100 | Repetições Monte Carlo por cenário |

Os cenários estão no dicionário `SCENARIOS`. As sementes das execuções são `1000 + i`, com `i` de 0 a 99.

## Reprodutibilidade

- As sementes são fixas, portanto **os erros de posição e as demais métricas de acurácia devem coincidir** entre execuções, com pequenas diferenças possíveis conforme a versão do NumPy.
- Os **tempos por ciclo** (`t_cycle_us`) dependem da máquina e mudam a cada execução. Foram medidos em computador de uso geral, sem otimização, e não representam desempenho em hardware embarcado.

## Limitações

- Simulação planar, com ruído gaussiano e uma única forma de trajetória.
- Na condição de referência, o EKF usa parâmetros de ruído coerentes com o gerador dos dados; o efeito de descalibração foi explorado apenas para o ruído do GNSS (`sensibilidade.py`).
- Os parâmetros dos sensores foram escolhidos pelos autores e não vêm de folhas de dados nem de medições em campo.
- Não há controlador de veículo: a camada de decisão da arquitetura é representada apenas pelo indicador de modo degradado.
- Os resultados indicam tendências e não permitem afirmar desempenho em campo nem economia de energia.

## Relação com a arquitetura do artigo

| Camada | Onde aparece no código |
|---|---|
| C1: adaptadores de sensores | Geração e sincronização das medições (dicionários `gps` e `whl`, vetores `wm` e `am`) |
| C2: estimador estocástico | Laço de predição e atualização do EKF em `run_once` |
| C3: monitor de consistência | Cálculo do NIS, validação, recuperação de travamento e r95 |
| C4: decisão e registro | Não implementada; representada por `pct_degradado` e `latencia_alerta_s` |

## Citação

Se usar este código, cite o artigo correspondente (referência completa a ser preenchida após a aceitação).
