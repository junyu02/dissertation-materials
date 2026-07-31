# Frozen stimulus library

`frozen_advice.csv` is the complete frozen advice library as deployed: 20
rows (10 fictional trader personas, each surfacing twice across the five
simulated months), exported directly from the platform's frozen database
table.

Columns: `trader` (fictional persona name), `month`, `round`, and the
deployed text fields — `recommendation_text` (the conclusion shown
byte-identically to both arms), `reasoning_text` (the primary frozen
reasoning body shown to the reasoning arm; mean 162 words), and the
direction- and holding-state-specific variants (`hold_text`,
`avoid_text`, `avoid_reasoning_text`, `unfollow_text`,
`unfollow_reasoning_text`).

All content is fabricated experimental stimulus material. The traders
are fictional personas, not real people and not participants; the market
events and returns are scripted. No participant data appears in this
directory.
