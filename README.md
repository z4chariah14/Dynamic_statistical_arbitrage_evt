# Dynamic Statistical Arbitrage with Tail Risk Control

I built this project to test whether a changing hedge ratio and an extreme value risk rule improve a relative value strategy between E mini S&P 500 and E mini Nasdaq 100 futures.

The main question is simple. If the relationship between ES and NQ changes over time, can a Kalman filter track that change better than one fixed hedge ratio? If it can, does an EVT risk rule protect the strategy when the estimated tail risk becomes unusually large?

## What the script does

Project.R downloads daily ES and NQ futures data, calculates log prices and returns, and checks the usual properties of financial time series. It tests for nonnormal returns, volatility clustering, unit roots, and cointegration.

The first strategy uses one Total Least Squares hedge ratio estimated during the formation period. The second strategy uses a Kalman filter to update the intercept and hedge ratio as new observations arrive.

The Kalman innovations are passed through either a standard GARCH model or a GJR GARCH model. The script chooses between them using BIC. It then fits a rolling Generalized Pareto Distribution to large absolute standardized innovations and estimates conditional Value at Risk and Expected Shortfall.

All trading decisions in the test period use information available at that point in time. Positions are applied one day later, and the backtest includes transaction costs.

## What I found

I used 2,203 daily observations from January 2018 through October 2, 2026. The first 1,321 observations were used for formation and the remaining 882 observations were reserved for testing.

ES and NQ returns were clearly nonnormal and showed strong volatility clustering. This supported using a conditional volatility model.

The static relationship was much weaker than I expected. The TLS spread failed its stationarity test, and both the Phillips Ouliaris and Johansen tests failed to find convincing evidence of cointegration. The static strategy entered one position and remained active for 495 test days because the spread did not reliably return to its formation period mean.

The Kalman innovations were stationary in both the formation and test periods. The dynamic strategy also had a much smaller drawdown than the static strategy. However, it entered only 13 positions and was active for 17 test days, so the small positive return is not enough evidence to claim a dependable trading advantage.

The EVT rule lowered annual volatility from 0.666 percent to 0.467 percent. It did not improve maximum drawdown or realized Expected Shortfall, and it changed the small positive return into a negative return.

The current results do not support the original idea that adding an EVT overlay automatically makes the strategy better. The stronger result is that a fixed ES and NQ relationship was unreliable, while dynamic hedging controlled the spread more effectively but produced too few trades to establish profitability.

## Out of sample results

The static TLS strategy returned negative 0.224 percent per year with a Sharpe ratio of negative 0.082 and a maximum drawdown of negative 5.392 percent.

The Kalman strategy returned 0.137 percent per year with a Sharpe ratio of 0.208 and a maximum drawdown of negative 0.853 percent.

The Kalman strategy with EVT risk control returned negative 0.137 percent per year with a Sharpe ratio of negative 0.291 and a maximum drawdown of negative 0.853 percent.

## Running the project

The script requires R and the following packages.

~~~r
install.packages(c(
  "quantmod",
  "tseries",
  "urca",
  "rugarch",
  "xts",
  "zoo"
))
~~~

Run the analysis from R with:

~~~r
source("Project.R")
~~~

The main numerical outputs are stored in project_results. The script also prints the model estimates and performance table and displays four diagnostic plots.

## Data limitation

The default data comes from Yahoo Finance. These daily continuous futures series are useful for a reproducible first study, but they are not properly roll adjusted contract data. The results should be treated as a model comparison rather than evidence of an executable futures strategy.

The script can also read a local CSV with Date, ES, and NQ columns by setting config$local_data_file.
