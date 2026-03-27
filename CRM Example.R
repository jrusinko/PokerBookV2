library(ggplot2)
library(reshape2)

# ============================================================
# USER PARAMETERS (EDIT ONLY THIS BLOCK)
# ============================================================

cards    <- 1:5
ante     <- 2
bet_size <- 2
n_iter   <- 10000

init_mode <- "random"   # "uniform" or "random"
seed      <- 123         # set NULL for no fixed seed

# ============================================================
# Setup
# ============================================================

if (!is.null(seed)) set.seed(seed)

n_cards <- length(cards)

pot_check <- 2 * ante
pot_bet   <- 2 * (ante + bet_size)
pot_fold  <- 2 * ante

# ============================================================
# Payoff Functions
# ============================================================

SD <- function(k, j, pot) {
  if (k > j) pot else -pot
}

# ============================================================
# Counterfactual Values (Player 1)
# ============================================================

v1_C <- function(i, tau) {
  k <- cards[i]
  opp_idx <- setdiff(seq_along(cards), i)

  mean(sapply(opp_idx, function(j_idx) {
    j <- cards[j_idx]
    (1 - tau[j_idx]) * SD(k, j, pot_check) +
      tau[j_idx] * (-pot_fold)
  }))
}

v1_B <- function(i, tau) {
  k <- cards[i]
  opp_idx <- setdiff(seq_along(cards), i)

  mean(sapply(opp_idx, function(j_idx) {
    j <- cards[j_idx]
    (1 - tau[j_idx]) * pot_fold +
      tau[j_idx] * SD(k, j, pot_bet)
  }))
}

# ============================================================
# Counterfactual Values (Player 2)
# ============================================================

v2_C <- function(j_idx, sigma) {
  j <- cards[j_idx]
  opp_idx <- setdiff(seq_along(cards), j_idx)

  mean(sapply(opp_idx, function(i_idx) {
    k <- cards[i_idx]
    showdown <- if (j > k) pot_check else -pot_check

    (1 - sigma[i_idx]) * showdown +
      sigma[i_idx] * (-pot_fold)
  }))
}

v2_B <- function(j_idx, sigma) {
  j <- cards[j_idx]
  opp_idx <- setdiff(seq_along(cards), j_idx)

  mean(sapply(opp_idx, function(i_idx) {
    k <- cards[i_idx]
    showdown <- if (j > k) pot_bet else -pot_bet

    (1 - sigma[i_idx]) * pot_fold +
      sigma[i_idx] * showdown
  }))
}

# ============================================================
# Regret Matching
# ============================================================

bet_prob <- function(RB, RC) {
  RB_pos <- max(RB, 0)
  RC_pos <- max(RC, 0)

  if (RB_pos + RC_pos > 0) {
    RB_pos / (RB_pos + RC_pos)
  } else {
    0.5
  }
}

# ============================================================
# Initialization
# ============================================================

if (init_mode == "uniform") {
  sigma <- rep(0.5, n_cards)
  tau   <- rep(0.5, n_cards)
} else if (init_mode == "random") {
  sigma <- runif(n_cards, 0, 1)
  tau   <- runif(n_cards, 0, 1)
} else {
  stop("init_mode must be 'uniform' or 'random'")
}

R1_B <- rep(0, n_cards)
R1_C <- rep(0, n_cards)

R2_B <- rep(0, n_cards)
R2_C <- rep(0, n_cards)

sigma_history     <- matrix(0, n_iter, n_cards)
sigma_avg_history <- matrix(0, n_iter, n_cards)

# ============================================================
# Main CFR Loop
# ============================================================

for (t in 1:n_iter) {

  sigma_history[t, ] <- sigma

  # --- Player 1 regrets ---
  for (i in seq_along(cards)) {

    vB <- v1_B(i, tau)
    vC <- v1_C(i, tau)

    v_current <- sigma[i] * vB + (1 - sigma[i]) * vC

    R1_B[i] <- R1_B[i] + (vB - v_current)
    R1_C[i] <- R1_C[i] + (vC - v_current)
  }

  # --- Player 2 regrets ---
  for (j in seq_along(cards)) {

    vB <- v2_B(j, sigma)
    vC <- v2_C(j, sigma)

    v_current <- tau[j] * vB + (1 - tau[j]) * vC

    R2_B[j] <- R2_B[j] + (vB - v_current)
    R2_C[j] <- R2_C[j] + (vC - v_current)
  }

  # --- Update strategies ---
  for (i in seq_along(cards)) {
    sigma[i] <- bet_prob(R1_B[i], R1_C[i])
  }

  for (j in seq_along(cards)) {
    tau[j] <- bet_prob(R2_B[j], R2_C[j])
  }

  # --- Running average ---
  sigma_avg_history[t, ] <- colMeans(sigma_history[1:t, , drop = FALSE])
}

# ============================================================
# Prepare Plot
# ============================================================

df <- as.data.frame(sigma_avg_history)
colnames(df) <- paste0("k=", cards)
df$Iteration <- 1:n_iter

df_long <- melt(df, id.vars = "Iteration",
                variable.name = "Card",
                value.name = "BetProb")

# ============================================================
# Plot
# ============================================================

ggplot(df_long, aes(x = Iteration, y = BetProb, color = Card)) +
  geom_line(linewidth = 1) +
  labs(
    title = paste("Average Betting Frequencies (init =", init_mode, ")"),
    x = "Iteration",
    y = "Probability of Betting",
    color = "Card"
  ) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal(base_size = 13)