---
  title: "Detecting Cheating and Collusion"
chapter-id: "detecting-cheating"
part: "Part IV: Information, Inference, and Deception"
---

  # Chapter Overview

  ## Stage 1: Chapter Scope

  **Chosen theme:** Statistical detection of cheating and collusion in poker using hypothesis testing and likelihood ratio methods.

**Central question:** Given observed data from a poker game, how can we determine — with quantifiable confidence — whether a player or group of players is behaving consistently with honest play, or whether their results are more likely explained by cheating or collusion?

  **Connection to poker:** Poker involves repeated decisions under uncertainty, making it a natural setting for statistical inference. Cheating and collusion distort the probability distributions that govern outcomes: a colluding team shares hole cards, a player using a marked deck gains information unavailable to honest players. These distortions leave statistical fingerprints that can be detected using formal hypothesis testing.

**Connection to earlier chapters:** This chapter draws on the probability distributions introduced in earlier chapters (hand frequencies, pot odds, equity calculations) as models of honest play. Deviations from those models form the basis of our statistical tests. The concept of expected value, used throughout the book, becomes here a reference point: a cheater's EV consistently exceeds what honest strategy can explain.

**Core mathematical tools:**

- Statistical hypothesis testing ($H_0$ vs. $H_1$, Type I and Type II error)
- The $p$-value as a measure of evidence against a null hypothesis
- The likelihood ratio and the Neyman--Pearson lemma
- Modeling poker outcomes as sequences of independent (or conditionally independent) observations
- Basic properties of the chi-squared and normal distributions (for large-sample approximations)

---

# Motivating Example

Suppose you are reviewing hand histories from an online poker platform. One player, call them Player X, has played 10,000 hands of no-limit hold'em over the past month. Over this sample, Player X has won money in an amount that seems extraordinary: their observed win rate is 12 big blinds per 100 hands (12 bb/100), roughly three times higher than the best documented win rates of elite professional players in comparable games.

Is this evidence of cheating? Maybe. But big samples can produce surprising results by chance, and a single summary statistic like win rate does not tell the whole story. To make a rigorous claim, we need a formal framework.

Now add another layer. Platform analysts notice that Player X frequently plays at the same tables as Player Y. When X and Y are both in a hand, X wins at an even higher rate. In hands where they are heads-up against each other, they exhibit an unusual pattern: one player almost always folds quickly, even when pot odds suggest calling is profitable. Chip transfers between the two accounts are suspiciously smooth.

This is the shape of a collusion problem. Two players coordinate — by sharing hole cards through an external channel, by soft-playing each other, or by signaling — to extract money from honest opponents. Their individual results may not look alarming in isolation, but their *joint* behavior is inconsistent with independent honest play.

Both scenarios — the individual's anomalous win rate and the pair's suspicious coordination — demand the same question: **How unlikely is this data under the assumption of honest play, and at what point does that unlikelihood constitute evidence of cheating?**

  The mathematics of hypothesis testing gives us a principled answer.

---

  # Mathematical Framework and Poker Theory

  ## Definitions and Setup

  We model poker outcomes as a probability space. Let each hand be an observation, and let the collection of hands played by one or more players constitute our dataset. Honest play defines a reference distribution; cheating or collusion distorts it.

### Definition 1: Statistical Hypothesis Test

**Definition.** A *statistical hypothesis test* is a procedure that uses observed data to evaluate a claim about an underlying probability distribution. The test specifies:

  1. A *null hypothesis* $H_0$, representing the baseline or default model (here, honest play);
2. An *alternative hypothesis* $H_1$, representing the competing model (cheating or collusion);
3. A *test statistic* $T$, a function of the data that measures departure from $H_0$;
4. A *rejection region* $\mathcal{R}$, a set of values of $T$ that lead us to reject $H_0$ in favor of $H_1$.

The test *rejects* $H_0$ when $T \in \mathcal{R}$, and *fails to reject* $H_0$ otherwise. Failing to reject $H_0$ does not establish that $H_0$ is true; it means only that the data do not provide sufficient evidence against it.

### Definition 2: Null and Alternative Hypotheses

**Definition.** In the context of cheating detection:

  - The *null hypothesis* $H_0$ is the claim that the player's outcomes arise from honest play — that is, from the probability distribution induced by the rules of the game and a legitimate strategy.
- The *alternative hypothesis* $H_1$ is the claim that the player's outcomes arise from a distribution inconsistent with honest play — for example, one where the player has access to hidden information.

**Type I error** (false positive): We reject $H_0$ when it is true — we accuse an honest player of cheating. The probability of this error is denoted $\alpha$, called the *significance level* of the test.

**Type II error** (false negative): We fail to reject $H_0$ when $H_1$ is true — a cheater goes undetected. The probability of this error is denoted $\beta$. The *power* of a test is $1 - \beta$.

In practice, platforms must choose $\alpha$ carefully. A very small $\alpha$ protects honest players from false accusations but allows more cheaters to escape detection. This is a genuine tradeoff, not a technicality.

### Definition 3: The $p$-value

**Definition.** The *$p$-value* of an observed test statistic $t$ is the probability, computed under $H_0$, of obtaining a test statistic at least as extreme as $t$:

  $$
  p\text{-value} = \mathbb{P}_{H_0}(T \geq t).
$$

  A small $p$-value means that the observed data would be very unlikely if $H_0$ were true. The $p$-value is *not* the probability that $H_0$ is true; it is the probability of seeing data this extreme or more extreme *given* that $H_0$ is true.

We reject $H_0$ when the $p$-value falls below the pre-specified significance level $\alpha$.

**Example.** If we model a player's per-hand profit under honest play as having mean $\mu_0$ and standard deviation $\sigma$, and we observe a sample mean $\bar{X}$ over $n$ hands, the test statistic

$$
Z = \frac{\bar{X} - \mu_0}{\sigma / \sqrt{n}}
$$

is approximately standard normal for large $n$. If the observed $z$-score is $z = 3.8$, the $p$-value for a one-sided test is approximately $\mathbb{P}(Z \geq 3.8) \approx 0.00007$. This is strong evidence against honest play.

### Definition 4: The Likelihood Ratio

**Definition.** Let $\mathbf{x} = (x_1, x_2, \ldots, x_n)$ denote the observed data. Suppose the data arise from a distribution with density (or probability mass function) $f(\mathbf{x} \mid \theta)$, where $\theta$ parameterizes the model. The *likelihood ratio* comparing hypothesis $H_1$ (parameter $\theta_1$) to hypothesis $H_0$ (parameter $\theta_0$) is

$$
\Lambda(\mathbf{x}) = \frac{f(\mathbf{x} \mid \theta_1)}{f(\mathbf{x} \mid \theta_0)}.
$$

A large value of $\Lambda(\mathbf{x})$ indicates that the data are much more probable under $H_1$ than under $H_0$, providing evidence in favor of $H_1$. The *log-likelihood ratio* is

$$
\ell(\mathbf{x}) = \log \Lambda(\mathbf{x}) = \log f(\mathbf{x} \mid \theta_1) - \log f(\mathbf{x} \mid \theta_0).
$$

---

## Main Mathematical Development

### Subsection A: Modeling Honest Play

Before we can detect deviation from honest play, we need a precise model of what honest play looks like.

Let $X_i$ denote the profit (in big blinds) of a player on hand $i$. Under honest play, the sequence $X_1, X_2, \ldots, X_n$ is assumed to be approximately i.i.d. with some distribution $F_0$, whose mean $\mu_0$ and variance $\sigma_0^2$ are determined by the player's strategy, the game structure, and the opponent pool.

In practice, $\mu_0$ is close to zero or slightly negative for most players (due to rake), and $\sigma_0^2$ is large due to the high variance of poker outcomes. Estimates of both can be derived from large databases of hand histories under verified honest play.

The key distributional fact for large $n$ is the **Central Limit Theorem**: regardless of the exact form of $F_0$,

$$
  \frac{\sqrt{n}(\bar{X} - \mu_0)}{\sigma_0} \xrightarrow{d} \mathcal{N}(0,1),
$$

  where $\bar{X} = \frac{1}{n} \sum_{i=1}^n X_i$. This justifies using the $Z$-test for large hand samples.

### Subsection B: The Neyman--Pearson Lemma and Optimal Detection

Among all tests with significance level at most $\alpha$, which test has the highest power? The **Neyman--Pearson Lemma** answers this question precisely.

**Theorem (Neyman--Pearson Lemma).** *Let $H_0 : \theta = \theta_0$ and $H_1 : \theta = \theta_1$ be simple hypotheses. The most powerful test of size $\alpha$ rejects $H_0$ when*

  $$
  \Lambda(\mathbf{x}) = \frac{f(\mathbf{x} \mid \theta_1)}{f(\mathbf{x} \mid \theta_0)} > k_\alpha,
$$

  *where $k_\alpha$ is chosen so that $\mathbb{P}_{H_0}(\Lambda(\mathbf{x}) > k_\alpha) = \alpha$.*

  **Proof sketch.** Consider any alternative test $\phi'$ also of size $\alpha$. Decompose the sample space into the region where $\phi'$ and the likelihood ratio test agree, and the regions where they disagree. On the region where the likelihood ratio test rejects but $\phi'$ does not, $\Lambda > k_\alpha$, so the likelihood ratio test gains more power per unit of type I error. A standard integration argument (or sum, in the discrete case) shows the likelihood ratio test achieves weakly higher power. $\square$

**Interpretation for cheating detection.** If we have a specific model for how a cheater behaves — say, a player who sees all opponents' hole cards has a known win-rate distribution $F_1$ — then the likelihood ratio test is the mathematically optimal way to distinguish cheaters from honest players. No other test of the same false-positive rate will catch more cheaters.

### Subsection C: Framing Detection as Distribution Discrimination

The central insight of this chapter can now be stated precisely.

**Result.** *Cheating detection can be framed as distinguishing between two competing probability distributions over observed game data. The likelihood ratio is the optimal statistic for this discrimination, in the sense of the Neyman--Pearson Lemma.*

  To apply this in practice, we model:

  - $H_0$: The data $\mathbf{x}$ arise from the distribution of honest play, $f(\mathbf{x} \mid \theta_0)$.
- $H_1$: The data $\mathbf{x}$ arise from the distribution of cheating, $f(\mathbf{x} \mid \theta_1)$.

For the collusion setting, suppose two players share their hole cards before acting. This grants them the ability to eliminate dominated strategies, fold losing hands they would otherwise call, and coordinate bet sizing. The effect on the distribution of outcomes is predictable and quantifiable:

  - Their *fold frequency* in hands where they would normally call increases.
- Their *win rate* in contested pots against honest players increases.
- The *correlation* between their net profits across hands deviates from zero.

Let $\mathbf{x}^{(A)}$ and $\mathbf{x}^{(B)}$ denote the hand-by-hand profits for players $A$ and $B$. Under honest play, their profits in hands where they are not in the same pot are independent. Under collusion, their profits are positively correlated (they help each other win) even when only one is in the pot.

We can test for collusion by examining the **sample correlation**

  $$
  \hat{\rho} = \frac{\sum_{i=1}^n (x_i^{(A)} - \bar{X}^{(A)})(x_i^{(B)} - \bar{X}^{(B)})}{\sqrt{\sum_i (x_i^{(A)} - \bar{X}^{(A)})^2 \cdot \sum_i (x_i^{(B)} - \bar{X}^{(B)})^2}}.
$$

  Under $H_0$ (independent honest play), $\hat{\rho} \approx 0$ for large $n$. A statistically significant positive $\hat{\rho}$ is evidence of coordination.

### Subsection D: Sequential Testing and Real-Time Monitoring

In practice, we do not wait until $n = 10{,}000$ hands before testing. Platforms monitor hands as they are played, updating their evidence continuously. This requires a different framework: **sequential hypothesis testing**.

The **Sequential Probability Ratio Test (SPRT)**, due to Wald, makes a decision at the first time $n$ at which the running likelihood ratio crosses a boundary:

  $$
  \ell_n = \sum_{i=1}^n \log \frac{f(x_i \mid \theta_1)}{f(x_i \mid \theta_0)}.
$$

  - If $\ell_n \geq B$, reject $H_0$ (flag as cheater).
- If $\ell_n \leq A$, accept $H_0$ (clear as honest, for now).
- Otherwise, continue collecting data.

The boundaries $A$ and $B$ are chosen to satisfy pre-specified Type I and Type II error rates $\alpha$ and $\beta$:

  $$
  A \approx \log\frac{\beta}{1 - \alpha}, \qquad B \approx \log\frac{1-\beta}{\alpha}.
$$

  The SPRT is optimal in the sense that, among all sequential tests with the same error rates, it requires the fewest observations in expectation to reach a decision — a property that matters enormously in a real-time platform context.

---

  ## Poker Connection

  ### What the model says in poker terms

  The honest-play null hypothesis encodes everything we know about how poker works: the rules, the deck composition, the betting structure, and the range of strategies available to players who act only on publicly visible information and their own private cards. Any player who exceeds the win rates achievable under this model is either extraordinarily lucky, exploiting an opponent pool in an extreme way, or cheating.

The likelihood ratio test tells us: given this specific model of cheating (e.g., hole-card sharing), how much more likely is this data under the cheating model than the honest model? If the answer is "vastly more likely," we have grounds to act.

### A concrete worked example

Suppose an honest player in a particular game has expected win rate $\mu_0 = 2$ bb/100 with standard deviation $\sigma_0 = 80$ bb/100 (a typical high-variance game). After $n = 5{,}000$ hands, we observe a win rate of $\bar{X} = 8$ bb/100.

The test statistic is:

  $$
  Z = \frac{(\bar{X} - \mu_0)}{\sigma_0 / \sqrt{n}} = \frac{(8 - 2)}{80 / \sqrt{5000}} = \frac{6}{1.131} \approx 5.30.
$$

  The $p$-value for a one-sided test is $\mathbb{P}(Z \geq 5.30) \approx 5.8 \times 10^{-8}$. This is extraordinarily small: if this player were honest, there is roughly a 1-in-17-million chance we would see a win rate this high. Under a significance level of $\alpha = 0.001$, we reject $H_0$ decisively.

### Assumptions and limitations

The model assumes that hands are approximately i.i.d., which is an idealization. In practice, player strategies are adaptive, opponent pools change over time, and stakes vary. The true distribution $F_0$ is unknown and must be estimated. These complications mean that statistical flags should be treated as *evidence* warranting investigation, not as proof of guilt.

Additionally, the $p$-value gives no information about *what kind* of cheating is occurring. A player running good over 5,000 hands might be lucky. A player sharing hole cards, using a real-time assistance tool, or colluding with an accomplice all produce different fingerprints, and a well-designed detection system uses multiple tests simultaneously.

---

  # Beyond Poker

  The framework developed in this chapter belongs to a much broader tradition of **anomaly detection** and **distribution testing** that appears across many domains.

## Finance: Detecting insider trading

A trader with access to nonpublic information — a form of insider information analogous to seeing opponents' hole cards — will consistently make profitable trades that cannot be explained by their stated strategy or public information. Regulators apply hypothesis testing and likelihood ratio methods to trading records to identify patterns inconsistent with honest market participation. The null hypothesis is that trades are consistent with information available to the public; the alternative is that they reflect a hidden information advantage.

## Sports integrity

Monitoring for match-fixing in professional sports follows exactly the same logic. The null distribution is derived from the performance of players in honest competition; outcomes inconsistent with that distribution — unusually timed errors, improbable patterns of play in matches with large betting line movements — trigger likelihood ratio-based alerts. The challenge of distinguishing cheating from extraordinary luck or skill is precisely the Type I / Type II error tradeoff described above.

## Network security and fraud detection

Online platforms of all kinds face the problem of distinguishing legitimate user behavior from bots, coordinated attacks, and fraudulent accounts. Here, the "honest play" null hypothesis is the distribution of behavior produced by real users (click rates, session lengths, action patterns), and the alternative hypothesis is the distribution produced by automated agents or coordinating groups. Sequential hypothesis testing methods, including generalizations of the SPRT, are standard tools in this setting.

## The broader lesson

In each application, the same mathematical structure appears:

1. Define a reference distribution for legitimate behavior.
2. Define one or more alternative distributions for anomalous behavior.
3. Use the likelihood ratio to measure how much more likely the observed data are under the alternative.
4. Choose thresholds that balance the cost of false accusations against the cost of missed detections.

Poker, with its well-defined rules and large volumes of logged hand data, is an unusually clean setting in which to study these ideas. But the principle — that cheating leaves a statistical signature, and that signature can be detected by comparing likelihoods — extends far beyond the card table.

---

# Homework Problems

1. **Conceptual problem.** Explain in your own words the difference between a $p$-value and the probability that $H_0$ is true. Why is this distinction important in the context of cheating accusations? Give an example of how confusing the two could lead to an unjust outcome.

2. **Computational problem.** A player has completed $n = 2{,}000$ hands with an observed win rate of $\bar{X} = 10$ bb/100. Assume that under honest play, the win rate distribution has mean $\mu_0 = 1$ bb/100 and standard deviation $\sigma_0 = 90$ bb/100. Compute the $Z$-statistic and the corresponding one-sided $p$-value. Would you reject $H_0$ at the $\alpha = 0.05$ level? At the $\alpha = 0.001$ level? Interpret your answer.

3. **Proof or derivation problem.** Let $X_1, \ldots, X_n$ be i.i.d. observations from a normal distribution with unknown mean $\mu$ and known variance $\sigma^2$. Show that the likelihood ratio test of $H_0 : \mu = \mu_0$ against $H_1 : \mu = \mu_1 > \mu_0$ is equivalent to rejecting $H_0$ when $\bar{X} > c$ for some threshold $c$ depending on $\alpha$. (This shows the $Z$-test is optimal in this setting, by the Neyman--Pearson Lemma.)

4. **Modeling or interpretation problem.** Two players, $A$ and $B$, are suspected of chip-dumping: deliberately losing chips to each other to transfer funds. In 200 heads-up pots between them, $A$ wins 142 and $B$ wins 58. Assume that under honest play, each heads-up pot is won by either player with equal probability.

   (a) Identify $H_0$ and $H_1$.

   (b) Compute the $p$-value using the normal approximation to the binomial.

   (c) What assumptions does this model make? Are they reasonable for heads-up poker?

5. **Optional challenge problem.** Derive the SPRT boundaries $A$ and $B$ as functions of $\alpha$ and $\beta$ under the approximation that the log-likelihood ratio $\ell_n$ crosses the boundaries exactly (ignoring overshoot). Then explain intuitively why the SPRT minimizes the expected sample size among all sequential tests with the same error rates. (A full proof is not required; a careful intuitive argument is sufficient.)

---

# Revision Checklist

- [x] The chapter has a clear title and central question.
- [x] The motivating example genuinely motivates the mathematics.
- [x] Definitions appear before use.
- [x] Notation is stable and consistent.
- [x] The chapter develops one main mathematical theme.
- [x] Poker is treated as a mathematical model, not as strategy advice.
- [x] The beyond-poker section broadens the chapter without changing its core mathematics.
- [x] The homework problems are varied and aligned with the chapter.
- [x] All AI-generated material has been checked and revised by the group.

---

# Submission Notes

## AI Collaboration Log

| Task | How AI Was Used | What the Group Should Verify |
|---|---|---|
| Definitions | Drafted formal definitions for hypothesis test, $p$-value, and likelihood ratio | Verify consistency with notation used elsewhere in the book |
| Neyman--Pearson Lemma | Generated theorem statement and proof sketch | Check proof sketch rigor; expand if a full proof is required |
| Worked example ($Z$-test) | AI calculated test statistic and $p$-value | Re-derive the $Z$-statistic computation and confirm the $p$-value using a table or software |
| Collusion model (correlation) | AI proposed correlation-based test | Verify the formula for $\hat{\rho}$; consider whether a more refined test is appropriate |
| SPRT section | AI drafted exposition and boundary formulas | Verify boundary approximations against Wald's original derivation |
  | Beyond-poker analogies | AI suggested finance, sports, and network security | Confirm that the analogies are accurate; add citations if needed |
  | Homework problems | AI generated candidate problems | Solve each problem independently to verify it is tractable with chapter content |
