# Changelog

## [1.8.0](https://github.com/Kohulan/ChemAudit/compare/v1.7.0...v1.8.0) (2026-06-18)


### Features

* **analytics:** make Butina clustering molecule cap deployment-configurable ([d84c3fd](https://github.com/Kohulan/ChemAudit/commit/d84c3fd614dd94e3aa4951b7644f0b20c4ee983a))
* **profiler:** enable SCScore via vendored weights + self-contained inference ([009d221](https://github.com/Kohulan/ChemAudit/commit/009d221e7a0bd3e67f89d1844952f7d4ad1e2ffe))
* **security:** add opt-in HMAC-signed admin tokens for replay resistance ([360ef54](https://github.com/Kohulan/ChemAudit/commit/360ef547e65b581f02d9b8843e8c86d7edcd5b5d))
* **theme:** add shared theme-aware chart color module ([4f25f2c](https://github.com/Kohulan/ChemAudit/commit/4f25f2c63dd27343a220a3e215f10199f59d929e))
* **ui:** replace native window.prompt/confirm with accessible, testable React modals ([0df42d0](https://github.com/Kohulan/ChemAudit/commit/0df42d01d013b6f2d16a1350883aa985c7d6f48f))


### Bug Fixes

* **a11y:** add dialog semantics and close-button label to ExportDialog ([0bb1bd2](https://github.com/Kohulan/ChemAudit/commit/0bb1bd2c59377d77df08f0f72f632964e621dd04))
* **a11y:** distinguish success/warning badges by icon and hue, not color alone ([ed66a07](https://github.com/Kohulan/ChemAudit/commit/ed66a073e430676ed1503ede138f33d5366d029b))
* **a11y:** keyboard-operable switch, modal semantics, tooltip associations ([497b703](https://github.com/Kohulan/ChemAudit/commit/497b703ffe1a9e4399e7242da0ab840bf618a7e6))
* **a11y:** raise text-muted contrast to WCAG AA in both themes ([41b7760](https://github.com/Kohulan/ChemAudit/commit/41b77603ddf0d6373b44e93690b1e01410e82dbe))
* **about:** address adversarial review of the distill pass ([646d48d](https://github.com/Kohulan/ChemAudit/commit/646d48da6c64cd1aa75222136205395f73bd6204))
* accessibility P0, palette stragglers, and responsive corners from re-audit ([e2cb165](https://github.com/Kohulan/ChemAudit/commit/e2cb165f9c2dddb18567dfc9d477525574ea0a11))
* address adversarial review findings from remediation pass ([c9c49a2](https://github.com/Kohulan/ChemAudit/commit/c9c49a2f5dbbb406a0a4d1f83588b91352257e75))
* **batch:** make ProgressTracker lazy Redis init thread-safe ([a06f400](https://github.com/Kohulan/ChemAudit/commit/a06f40083ae35a5003114d642c400076047b8672))
* **batch:** narrow Celery autoretry_for to transient errors only ([693ff3b](https://github.com/Kohulan/ChemAudit/commit/693ff3b7118800161893d40326b53b82584e2046))
* **celery:** add global task time limits to prevent worker starvation ([5bce2c1](https://github.com/Kohulan/ChemAudit/commit/5bce2c122473496b0a45279d55ceb0b88c28423a))
* **deps:** bump dompurify &gt;=3.4.11 and patch frontend audit findings ([c86bc0f](https://github.com/Kohulan/ChemAudit/commit/c86bc0fd84ad301c5091201379823718a9813578))
* **deps:** bump dompurify &gt;=3.4.11 and patch frontend audit findings ([a16142c](https://github.com/Kohulan/ChemAudit/commit/a16142c1234db4b511d7c1cc6940e0dadb814ac4))
* **deps:** bump joi to &gt;=18.2.1 in docs-site (CVE-2026-48038) ([5324d8e](https://github.com/Kohulan/ChemAudit/commit/5324d8ee186d211e83e1258e0aad0a1803e81c6e))
* **deps:** bump joi to &gt;=18.2.1 in docs-site (CVE-2026-48038) ([e4e7a4a](https://github.com/Kohulan/ChemAudit/commit/e4e7a4ac84a1c07f2332e6caa158513a256f61d3))
* **deps:** cap fastapi &lt;0.137 to avoid instrumentator route crash ([053ff43](https://github.com/Kohulan/ChemAudit/commit/053ff436f6e12c22cca33edb5059d87ab62a739f))
* **deps:** cap fastapi &lt;0.137 to avoid instrumentator route crash ([937b7f0](https://github.com/Kohulan/ChemAudit/commit/937b7f08f9b62a0515cca664267e79892601b513))
* **deps:** force esbuild &gt;=0.28.1 in frontend (GHSA-gv7w-rqvm-qjhr, GHSA-g7r4-m6w7-qqqr) ([35c8b41](https://github.com/Kohulan/ChemAudit/commit/35c8b41c7a79ff592de2a391a261382d3dbcc3c2))
* **deps:** force esbuild &gt;=0.28.1 in frontend (GHSA-gv7w-rqvm-qjhr, GHSA-g7r4-m6w7-qqqr) ([d78ce28](https://github.com/Kohulan/ChemAudit/commit/d78ce28a7c5ec2808dfb2d8adb8433e0da5aafd8))
* **deps:** fully resolve js-yaml alert by moving gray-matter to js-yaml 4 ([0a347c1](https://github.com/Kohulan/ChemAudit/commit/0a347c1eaf3d131bc37c25e8b8d0bdbd05ec58bf))
* **deps:** patch docs-site Dependabot alerts via npm overrides ([f74abf4](https://github.com/Kohulan/ChemAudit/commit/f74abf4723d83e55c42e64c477aef0f436d859c2))
* **deps:** resolve Dependabot security alerts in frontend and docs-site ([f81646a](https://github.com/Kohulan/ChemAudit/commit/f81646ac4b59cf5f911fb99b903060d93ffe6782))
* **deps:** resolve Dependabot security alerts in frontend and docs-site ([236dfbd](https://github.com/Kohulan/ChemAudit/commit/236dfbdc295abd1255a6363bf0ed3dfc06c14111))
* **deps:** sync frontend lockfile after vitest 4 bump ([00f17c0](https://github.com/Kohulan/ChemAudit/commit/00f17c0902e03026a63cebb0cbcc727e2566fd86))
* **deps:** sync frontend lockfile after vitest 4 bump ([9c155dc](https://github.com/Kohulan/ChemAudit/commit/9c155dc743850027c66fcc269529367d9c39d7ad))
* frontend audit remediation (a11y, performance, responsive, theming) ([e1e71ff](https://github.com/Kohulan/ChemAudit/commit/e1e71ff35ad626729f7d9a95ad1eba45c421841d))
* **iupac:** skip OPSIN re-init after failure to avoid inconsistent JVM state ([03b6708](https://github.com/Kohulan/ChemAudit/commit/03b670867d5719b815c167c5e0c5f91045718a33))
* **profiler:** make SCScore unavailability observable + document NumPy 2.x enablement ([11abdc6](https://github.com/Kohulan/ChemAudit/commit/11abdc6daa278a60279d0e501ed11418b805742a))
* **profiler:** validate SYBA SMILES input and make subprocess shell=False explicit ([982e9cd](https://github.com/Kohulan/ChemAudit/commit/982e9cd4e7d90f3ea1c74abd3942696e56978aac))
* **redis:** bound memory and use volatile-lru so broker messages are never evicted ([1f0ec2e](https://github.com/Kohulan/ChemAudit/commit/1f0ec2ea311c0a9ca5b37be2048c2a6d3228ea83))
* regressions found by adversarial sweep of the cleanup branch ([bdf0b75](https://github.com/Kohulan/ChemAudit/commit/bdf0b7525ec5527b8e7e8a0033278300af8953b5))
* **responsive:** add mobile breakpoints to stat and comparison grids ([daf3d0b](https://github.com/Kohulan/ChemAudit/commit/daf3d0b2d901667bddd1640e7fcaeb457afde572))
* **responsive:** enforce 44px touch target on mobile menu button ([0ea6ed8](https://github.com/Kohulan/ChemAudit/commit/0ea6ed8dcc3896e08d0b9f71145fdb0041ef8519))
* **responsive:** mobile-safe dropdown and SMILES column, intermediate grid breakpoints ([b84d09e](https://github.com/Kohulan/ChemAudit/commit/b84d09e9551edf41bea2cc280e1e3e1300b936e1))
* **scoring:** guard module-level sascorer imports against missing RDKit Contrib ([e1afd70](https://github.com/Kohulan/ChemAudit/commit/e1afd70bac5f9db2a36c48771343536944293637))
* **security:** escalate insecure-default warnings and add ALLOW_INSECURE_DEFAULTS gate ([7f4a1cf](https://github.com/Kohulan/ChemAudit/commit/7f4a1cf48f14344ff92e99caec0ed091d2fd11d9))
* **security:** stop logging secret-named settings; document HMAC-SHA256 intent ([b547823](https://github.com/Kohulan/ChemAudit/commit/b547823dfc090d046280b2a9586d917417b35ddc))
* **theme:** AlertGroupCard severity colors use status tokens with dark variants ([ded142c](https://github.com/Kohulan/ChemAudit/commit/ded142c2eeff3a6fc0e6c063af67614d7f0f2544))
* **theme:** batch and shared components use tokens; success is amber not green ([b449847](https://github.com/Kohulan/ChemAudit/commit/b449847a1e47e0ff546d06aac2afad98d2170e22))
* **theme:** batch score histogram uses warm brand scale, not green ([df65c6e](https://github.com/Kohulan/ChemAudit/commit/df65c6ea20671b2691d38c763e3b7f37ebe40e8e))
* **theme:** integrations panels use surface/border tokens ([e4c2fe0](https://github.com/Kohulan/ChemAudit/commit/e4c2fe00b939cb0e7d327d2cc9951ca272b962bc))
* **theme:** legible molecule renderings in dark mode via inversion filter ([bb17acc](https://github.com/Kohulan/ChemAudit/commit/bb17acc5d089ba1d5cd258ece260dae8f0e3d7d7))
* **theme:** NP-likeness reads red (synthetic) to green (natural) ([ad44bb7](https://github.com/Kohulan/ChemAudit/commit/ad44bb7fffc33c0b29e857dd70afa2561946f15b))
* **theme:** replace decorative purple/indigo with brand palette (database identity colors kept) ([1a66ea8](https://github.com/Kohulan/ChemAudit/commit/1a66ea8801884acc401d958a19a69a0caaceff6f))
* **theme:** retire green-for-success and cool score palettes; badge contrast to AA ([2d2148b](https://github.com/Kohulan/ChemAudit/commit/2d2148bc7364d75622675c6f29fdd236360a936d))
* **theme:** standardization components use surface/border tokens ([d48c311](https://github.com/Kohulan/ChemAudit/commit/d48c311dc0ea0a30c6cd0f63386f365ab884a267))
* **theme:** warm gauge tracks and on-palette chart reference colors ([cb62e41](https://github.com/Kohulan/ChemAudit/commit/cb62e417de81077a081b0bf5cde74a9761e8db09))
* **ui:** anchor auto-compare knob to the track's left edge ([f7028ec](https://github.com/Kohulan/ChemAudit/commit/f7028ecc92c9478073a57305081ca5efa7dbbd0f))
* **ui:** inline delete confirmation, lucide icon swaps, copy and asset cleanup ([01ae4c8](https://github.com/Kohulan/ChemAudit/commit/01ae4c8ae0d20491211dd204227f3464d8bff63e))
* **ux:** scroll to top on link navigation app-wide ([a7063c3](https://github.com/Kohulan/ChemAudit/commit/a7063c3875cd9e687a6db6ffc7eb22d4afe7f8aa))
* **websockets:** prevent duplicate broadcasts from superseded subscriber + add WS test suite ([4cbfb5e](https://github.com/Kohulan/ChemAudit/commit/4cbfb5ebb2448eccb41c67cbd7d2daadae5d2c06))
* **websockets:** use redis aclose() instead of deprecated close() ([b6ca465](https://github.com/Kohulan/ChemAudit/commit/b6ca4656f3b99b0259de47efca6fa016af0140f7))


### Performance Improvements

* **analytics:** isolate expensive analytics on dedicated queue + worker with time limits ([dc82cd1](https://github.com/Kohulan/ChemAudit/commit/dc82cd18c89ea7cec56deba7ede1a7433244fedb))
* **batch:** store results as Redis list for range-query pagination without full deserialize ([f75fcc0](https://github.com/Kohulan/ChemAudit/commit/f75fcc0ffe770cdb63b1952ad57d39d56f9b641f))
* **clustering:** vectorize Tanimoto distance matrix with BulkTanimotoSimilarity ([b6107bc](https://github.com/Kohulan/ChemAudit/commit/b6107bc6ad86ae1dde634a4270ac8b29d2a370f0))
* lazy-load scoring and standardization result panels on home page ([c7bb63d](https://github.com/Kohulan/ChemAudit/commit/c7bb63d751daf1479279241d1b136dd6aa3c6bea))
* opacity-only ambient orbs, capped blur radii, reduced-motion support ([ad5d20c](https://github.com/Kohulan/ChemAudit/commit/ad5d20cad8d1593cb20695e56de4ff83abdbc7cb))
* **profiler:** replace per-call SYBA subprocess with persistent worker (JSON pipe, restart, timeout) ([24057f5](https://github.com/Kohulan/ChemAudit/commit/24057f502bd99084601c8bfb886be7a4a322930a))
* **redis:** use pooled singleton clients instead of per-call connections ([30f0898](https://github.com/Kohulan/ChemAudit/commit/30f0898e5bc8a435c1187a3eaee6adac9f5f9928))
* route-gate RDKit WASM init; skip splash on structure-free pages ([ec20890](https://github.com/Kohulan/ChemAudit/commit/ec20890e454839e8a6b4094f46492f7ba4e30e85))
* **security:** pin RDKit loader version, add SRI, defer parsing ([29d627b](https://github.com/Kohulan/ChemAudit/commit/29d627be09807278fc814ef8e56c2ac00ba33ae5))
* serve resized logo variants instead of 764KB original ([8be0f74](https://github.com/Kohulan/ChemAudit/commit/8be0f74a85586806c0c6b2799d71f9cdd4f590d3))
* skip offscreen row rendering in batch results table ([3ead8ec](https://github.com/Kohulan/ChemAudit/commit/3ead8ec5490d7f48ab170d4825578079f350aed2))
* split recharts, framer-motion, gsap, and react vendor chunks ([b0c9a10](https://github.com/Kohulan/ChemAudit/commit/b0c9a10f201eba3c0347e5583b38403525e7f672))
* **validation:** run blocking Celery task.get in executor to free event loop ([a7a082a](https://github.com/Kohulan/ChemAudit/commit/a7a082a2bc32046be9c31883e963d45d577f5943))


### Documentation

* add frontend audit remediation plan ([224fa7d](https://github.com/Kohulan/ChemAudit/commit/224fa7d6438223c2587528c0c42cc4545f766ca8))
* **plan:** mark 2.5 SCScore done; all audit items now resolved ([c790061](https://github.com/Kohulan/ChemAudit/commit/c790061c5680d0f4b2ca69a2def72e64ae05f5d0))
* **plan:** mark 4.1 complete; reconcile 2.6/3.3/3.5 as done (were stale 'deferred') ([e611777](https://github.com/Kohulan/ChemAudit/commit/e611777efe102fdad5bc974e9bd0e04a29449588))
* **plan:** record [#7](https://github.com/Kohulan/ChemAudit/issues/7) documented + [#21](https://github.com/Kohulan/ChemAudit/issues/21) made configurable (gaps closed) ([7affa7d](https://github.com/Kohulan/ChemAudit/commit/7affa7d9baf100189d16aabaa56794e3216e56be))
* **plan:** record Wave 3 + Wave 4.3 completion and remaining-item status ([2654fa0](https://github.com/Kohulan/ChemAudit/commit/2654fa08f509dcb8001667eaf1418a40e9831c8d))
* **plan:** record Wave 4.2/4.4/4.5 completion + real-stack validation ([c88f583](https://github.com/Kohulan/ChemAudit/commit/c88f583e8453e96f628da4bcb8e3022853837793))
* **plan:** record WS fail-open decision and Wave 2 completion status ([91e6e1f](https://github.com/Kohulan/ChemAudit/commit/91e6e1f3aa8dd80fa6229ee9392314751cb3ffdd))
* **ws:** document deliberate fail-open WebSocket ownership posture ([b644e5a](https://github.com/Kohulan/ChemAudit/commit/b644e5af5801414abb1024d87d3cf5521ae02236))

## [1.7.0](https://github.com/Kohulan/ChemAudit/compare/v1.6.0...v1.7.0) (2026-05-11)


### Features

* **delight:** bench-register batch completion language (impeccable D3) ([5e061a7](https://github.com/Kohulan/ChemAudit/commit/5e061a78535d9415f3d35c1dbdb61d0797f57164))
* **delight:** chemistry-accurate loading messages cycle (impeccable D1) ([482c46d](https://github.com/Kohulan/ChemAudit/commit/482c46d0a1d00154409a7c756124f2f88ec2bf9b))
* **delight:** chemistry-aware empty state for structure preview (impeccable D2) ([97a95ef](https://github.com/Kohulan/ChemAudit/commit/97a95ef8a56e3c33f35a8d60d68fc0d6dc973b02))
* **delight:** console citation easter egg for developer-chemists (impeccable A4) ([b4ae33a](https://github.com/Kohulan/ChemAudit/commit/b4ae33a4fd6ff425b59abe54235c033f775f64fb))
* **delight:** structure preview caption shows formula and MW (impeccable D4) ([6302bbf](https://github.com/Kohulan/ChemAudit/commit/6302bbf62df3334578192513a97292548c459f72))
* **integrations:** add QLever as fallback Wikidata SPARQL endpoint ([4103627](https://github.com/Kohulan/ChemAudit/commit/41036276037c300902806d93226954773ffc02df))
* **ui:** All Checks expands to full-width 4-up grid ([19f23f2](https://github.com/Kohulan/ChemAudit/commit/19f23f2b214843a199b962be221e5d3e2f4fecc4))
* **ui:** All Checks panel morphs from right column to full-width on expand ([35747b8](https://github.com/Kohulan/ChemAudit/commit/35747b8b0503a4deae8381529df4c2a1ecc456b3))
* **ui:** categorize All Checks panel into 7 sections with hover delight ([70b448a](https://github.com/Kohulan/ChemAudit/commit/70b448a23186009679e59cd77340be61cf6d3885))
* **ui:** colorize Molecule Info stat tiles + input-type badge (impeccable colorize) ([9c21fe0](https://github.com/Kohulan/ChemAudit/commit/9c21fe069b888af4c18fb90635e3c8d0fefcb5c6))
* **ui:** diagnostic tab dot colours by result state (impeccable E1) ([90dfc2c](https://github.com/Kohulan/ChemAudit/commit/90dfc2c1750066b0d301144ff3ac120c95d669aa))
* **ui:** distill SingleValidation to 3 primary tabs + disclosure (impeccable phase 4) ([6254bce](https://github.com/Kohulan/ChemAudit/commit/6254bced8f5e6ac5b10d458fdca73a56a9bb4739))
* **ui:** parse error sub-typing with bench-grade hints (impeccable E2) ([e32ec0f](https://github.com/Kohulan/ChemAudit/commit/e32ec0fdf02c5b2a2608215ed35387ea0c6c26f6))
* **ui:** semantic icon colour on action row buttons (impeccable colorize) ([d73c8dc](https://github.com/Kohulan/ChemAudit/commit/d73c8dc71689ce85f142e37b29c880fd666aa661))
* **ui:** single floating All Checks card with sequenced position+size animation ([30604c3](https://github.com/Kohulan/ChemAudit/commit/30604c3ac5c95be336ee0e157f119273303f0481))
* **ui:** tab has-result dot indicator (impeccable critique P1) ([28cdf0f](https://github.com/Kohulan/ChemAudit/commit/28cdf0fe9da7f0777115f1942befd85ce4c0dfe8))


### Bug Fixes

* **a11y:** hasResult dot announces 'has results' to screen readers (impeccable A3) ([d9c989a](https://github.com/Kohulan/ChemAudit/commit/d9c989a811d430d6a1e694791802d82292c97020))
* **deps:** resolve 15 dependency vulnerabilities across frontend and docs-site ([3cc5321](https://github.com/Kohulan/ChemAudit/commit/3cc53213255d0441b1866f2c455c8d76d706087f))
* **lint:** remove unused no-console eslint-disable directive ([46ae427](https://github.com/Kohulan/ChemAudit/commit/46ae42797dbc5de3fbc4258592ec1c6e0217a04e))
* **ui:** All Checks card moves straight down (not diagonally) ([5ce0536](https://github.com/Kohulan/ChemAudit/commit/5ce053610e3a8055f3c3aef1de983de29460196d))
* **ui:** clarify error messages, empty states, and form labels ([2651d04](https://github.com/Kohulan/ChemAudit/commit/2651d047b36c3e0352a3c8fe8ab570042a2db38a))
* **ui:** clarify NP-Likeness copy (impeccable clarify NP-Likeness) ([e5352cc](https://github.com/Kohulan/ChemAudit/commit/e5352cc5e595a7aea94fd385d2cd63f07d74707a))
* **ui:** close 24px gap before All Checks and equalize score tile heights ([4a88c47](https://github.com/Kohulan/ChemAudit/commit/4a88c471c63cfe0a3215254337544fa48aba7224))
* **ui:** database lookup error now surfaces to user (impeccable B1, P0) ([49e1ef0](https://github.com/Kohulan/ChemAudit/commit/49e1ef005707bf5956dfc0b399327b4be25083d1))
* **ui:** error recovery + minor observations (impeccable phase 6) ([df7d375](https://github.com/Kohulan/ChemAudit/commit/df7d375a4dedda8e0c5deae5645f0df143efdc0c))
* **ui:** gate Safety / Compound Profile auto-screen on validation result (impeccable B2, P1) ([144eed0](https://github.com/Kohulan/ChemAudit/commit/144eed0ab0c90e18d44bca674d5e05d1ab243d59))
* **ui:** gate Score button on validation completion (impeccable critique P1) ([52d69bd](https://github.com/Kohulan/ChemAudit/commit/52d69bd4faa3b0faf5f169872c47c20f8a2dcba3))
* **ui:** gate validate-tab right-column panels on activeTab (impeccable B3, P1) ([c79be9e](https://github.com/Kohulan/ChemAudit/commit/c79be9e744f6f3e42c53ea886f902629f9891519))
* **ui:** inline check descriptions + bigger cell height in All Checks grid ([0d3cf4c](https://github.com/Kohulan/ChemAudit/commit/0d3cf4c8be2976bb43d51f50757a00f6976e7018))
* **ui:** IssueSeverityTags 'All clear' uses gold-amber per Warm-Status Rule (impeccable A2) ([af112a2](https://github.com/Kohulan/ChemAudit/commit/af112a2c290204d02868fbdaaab387e0f244d5ae))
* **ui:** make More analysis disclosure always discoverable on SingleValidation ([5a94186](https://github.com/Kohulan/ChemAudit/commit/5a941868430c6358fec9f4fe5485856bbff319ae))
* **ui:** measure All Checks bounds via offset chain ([3341823](https://github.com/Kohulan/ChemAudit/commit/33418239c57bffd45fd8362551014118acb1eb75))
* **ui:** measure All Checks bounds via offset chain to avoid transform jitter ([29b91b8](https://github.com/Kohulan/ChemAudit/commit/29b91b884fa2f7a779521a6469d684f4dc920900))
* **ui:** merge duplicate lucide-react import in SingleValidation ([94660d5](https://github.com/Kohulan/ChemAudit/commit/94660d5804e344ebc90ab03bc04a97ac8152708d))
* **ui:** move ScoringResults + DatabaseComparisonPanel into right column (impeccable C1, P1) ([3bf26ba](https://github.com/Kohulan/ChemAudit/commit/3bf26ba4e5c62059781b022535fe4957b5348288))
* **ui:** proper claymorphism on NP-Likeness track + marker (impeccable clarify NP-Likeness) ([0c148b4](https://github.com/Kohulan/ChemAudit/commit/0c148b495085c43ce0041055abfcd589ccb9e15c))
* **ui:** purge off-brand colors and absolute-ban anti-patterns ([00511e6](https://github.com/Kohulan/ChemAudit/commit/00511e6b4e0f546b40bea87ea491eb7127692a2d))
* **ui:** quiet About page background motion (impeccable phase 5) ([0e00eec](https://github.com/Kohulan/ChemAudit/commit/0e00eec271e266c90b1455855c12d633a280d20b))
* **ui:** replace amber-on-amber pass styling with prominent gold + CheckCircle2 (impeccable critique P0) ([fa6b36a](https://github.com/Kohulan/ChemAudit/commit/fa6b36a83fd0c89190d0fe2d18f8c87af911f533))
* **ui:** replace lock emoji with Lock icon (impeccable A1) ([3537002](https://github.com/Kohulan/ChemAudit/commit/35370025b5c3e2f1f8e64ccbe6788bbf9e4e4a75))
* **ui:** replace More analysis disclosure with two-row tab layout ([a6c2714](https://github.com/Kohulan/ChemAudit/commit/a6c27145470087752732f10738f1795dd4d86655))
* **ui:** replace tab descriptions with first-actionable guidance (impeccable F.P2-A) ([8cbe4ff](https://github.com/Kohulan/ChemAudit/commit/8cbe4ffae7cfe4564b2926feb8c4abf822ca192c))
* **ui:** residual color and absolute-ban cleanup (phase 2.5) ([ed63f02](https://github.com/Kohulan/ChemAudit/commit/ed63f0240a5040fac79353628e01f4c6a38fdc6b))
* **ui:** revert ScoringResults to full-width below grid + auto-scroll on land ([2c2b2bf](https://github.com/Kohulan/ChemAudit/commit/2c2b2bf468aa90d86763bce4c8af7699fc620597))
* **ui:** SingleValidation error card recovery parity (impeccable critique P0) ([17025ef](https://github.com/Kohulan/ChemAudit/commit/17025ef80aa6c6beda10e58da26b017d8bf72566))
* **ui:** smooth single-element morph for All Checks panel ([4b59703](https://github.com/Kohulan/ChemAudit/commit/4b59703ae1b2ad99a59e38e6302903601a00e931))
* **ui:** unify SingleValidation action row buttons (impeccable clarify) ([272c465](https://github.com/Kohulan/ChemAudit/commit/272c46561d01883203fba2c44e59a8b95daac255))
* **ui:** unify two SingleValidation tab rows into a single tray ([5605acf](https://github.com/Kohulan/ChemAudit/commit/5605acf93ac46f03dc633e8f10a2118635155c99))


### Documentation

* **ui:** clarify ScoreChart 'cool' variant is rose/crimson not blue (impeccable F.P2-C) ([23d8548](https://github.com/Kohulan/ChemAudit/commit/23d8548ada41599c49d093b6b8a52105a5f95b68))

## [1.6.0](https://github.com/Kohulan/ChemAudit/compare/v1.5.1...v1.6.0) (2026-04-24)


### Features

* **security:** log suspicious patterns in uploaded batch files ([fe72c6b](https://github.com/Kohulan/ChemAudit/commit/fe72c6b6b1a7e3b05a4756c13a49b3790d7942a3))
* **skills:** add ChemAudit Claude Code plugin with 8 agent skills ([af91d6a](https://github.com/Kohulan/ChemAudit/commit/af91d6a8617ea790cfd9ee0f48fcd178b1d45316))
* **skills:** add ChemAudit Claude Code plugin with 8 agent skills ([8dd8eb2](https://github.com/Kohulan/ChemAudit/commit/8dd8eb22f0c0da52a3ae87ef57a9ce35e582fd8a))
* **ui:** replace BB-8 toggle with chemistry-themed flask toggle ([9984a57](https://github.com/Kohulan/ChemAudit/commit/9984a5727ced739b42535ae4752033750c7d3038))


### Bug Fixes

* **analytics:** remove unreachable salt-form dedup fast-path ([914ab97](https://github.com/Kohulan/ChemAudit/commit/914ab973d9a0a2a6136821fab10e5d606b026c5d))
* **batch:** accept single-column and headerless delimited files ([c9d0033](https://github.com/Kohulan/ChemAudit/commit/c9d00331a3e78bef80c86b89569aad75cd4cee21))
* **batch:** accept single-column and headerless delimited files ([025f9bc](https://github.com/Kohulan/ChemAudit/commit/025f9bc31adacfbb8d232a358644a4ee1e714e28))
* **deploy:** make deploy.sh succeed on fresh installs ([94bb599](https://github.com/Kohulan/ChemAudit/commit/94bb5995074af7aa5daefad48a6f0ce1c45942e1))
* **deps:** resolve uuid and follow-redirects vulnerabilities in docs-site ([cf8a43d](https://github.com/Kohulan/ChemAudit/commit/cf8a43dbf3a63cfa82f08754820234d416494a84))
* **diagnostics:** highlight MCS-based structural diff and fix thumbnail overflow in cross-pipeline panel ([00c8aa2](https://github.com/Kohulan/ChemAudit/commit/00c8aa2029f209faf5c21f1687d65ecd2fb62458))
* **security:** run suspicious-content scan before strict validation ([f62a585](https://github.com/Kohulan/ChemAudit/commit/f62a585dd3b5cf28516424bb81ef5cfe207b9fa3))
* **security:** sanitize scscore error returned to client ([df13844](https://github.com/Kohulan/ChemAudit/commit/df13844b309337c13d7b7c67483d056315e7a838))
* **ui:** add info text above validate buttons in scoring and profile tabs ([66d72b9](https://github.com/Kohulan/ChemAudit/commit/66d72b9f45f6d529c4d3b7b18541c743d562d8d5))
* **ui:** prevent scoring profiles auto-run and add validate buttons ([d049375](https://github.com/Kohulan/ChemAudit/commit/d0493759f187946c39742172ef2521c95a6db250))
* **validation:** skip name resolver when input is valid SMILES ([0388336](https://github.com/Kohulan/ChemAudit/commit/038833681af2140583a0c3176b383ccb89bccf74))


### Documentation

* **api:** mark refreshCsrfToken as intentional public escape hatch ([2906e60](https://github.com/Kohulan/ChemAudit/commit/2906e60611d95273baa14b66a4eb80be77ecdcaf))

## [1.6.0](https://github.com/Kohulan/ChemAudit/compare/v1.5.0...v1.6.0) (2026-04-14)


### Features

* **analytics:** wire taxonomy bar chart click-to-filter with molecule indices ([c9231ee](https://github.com/Kohulan/ChemAudit/commit/c9231ee556f0d81765701ac1a0a5d7e61f57ca93))
* **diagnostics:** add duplicate column header detection to CSV pre-validator ([97b82ef](https://github.com/Kohulan/ChemAudit/commit/97b82efd04a32a425f7818dfd998b85753addb0e))
* **export:** implement PDF audit section rendering in batch report template ([3da093e](https://github.com/Kohulan/ChemAudit/commit/3da093e66f8d5759a67fdcdf6229005df228d46d))
* **integrations:** add SureChEMBL logo and About page acknowledgment ([613b60f](https://github.com/Kohulan/ChemAudit/commit/613b60fdf883e71b2b392b5b182c87b0cc75a501))
* **integrations:** add SureChEMBL to database lookup description list ([6c3e26c](https://github.com/Kohulan/ChemAudit/commit/6c3e26c6bd9f71d588aaa8bfbbd1ea46acdc3e43))
* **integrations:** wire SureChEMBL patent lookup into database results UI ([6bd4a46](https://github.com/Kohulan/ChemAudit/commit/6bd4a46ce07eeaf7639b5e084ae42b6421907b08))


### Bug Fixes

* **integrations:** update SureChEMBL source fallback and type comment ([3fd5862](https://github.com/Kohulan/ChemAudit/commit/3fd586205069b4a26e17270e6e6ced4d0a3be5ff))
* **qsar:** use correct prevalidation gate key in batch upload ([3d446e5](https://github.com/Kohulan/ChemAudit/commit/3d446e518b87d5ade1e69fe058cca5178e91c04d))
* sticky header and code cleanup ([8a28371](https://github.com/Kohulan/ChemAudit/commit/8a283717a9cbf5f9d0d95c5f3c62e27f981942c2))
* sticky header and code cleanup ([ae66b3f](https://github.com/Kohulan/ChemAudit/commit/ae66b3f2392bd8a2df5173679acba73797c95900))


### Documentation

* **mcp:** add MCP server documentation and fix nginx proxy for /mcp ([f2347b1](https://github.com/Kohulan/ChemAudit/commit/f2347b1f0ece5c0977ac43ecc87a182d843a08d1))
* **mcp:** add MCP server documentation and fix nginx proxy for /mcp ([977e386](https://github.com/Kohulan/ChemAudit/commit/977e38664072c4e94b5f1d5439d5a445c62e8832))
* rewrite README for v1.5.0 with complete feature coverage ([3575bda](https://github.com/Kohulan/ChemAudit/commit/3575bdae9820a62f5e9c3e10bba655eec27365af))

## [1.5.0](https://github.com/Kohulan/ChemAudit/compare/v1.4.0...v1.5.0) (2026-04-13)


### Features

* **about:** add developer profile card with claymorphism design ([861104a](https://github.com/Kohulan/ChemAudit/commit/861104a39e6617eccd02165141a489a48bc7937b))
* add Wikidata lookup, BB-8 theme toggle, clustering/taxonomy info boxes, update About page ([7ea2957](https://github.com/Kohulan/ChemAudit/commit/7ea2957c20feb5c51c1d9767ca23866d4893d6e7))
* **batch:** move Export All to top bar, Export Selected to results header ([194da7e](https://github.com/Kohulan/ChemAudit/commit/194da7e046a36e4022f9ea331aa0eec2b831ab10))
* **batch:** widen results table, improve info icons, enrich collision groups ([4e1b912](https://github.com/Kohulan/ChemAudit/commit/4e1b912e98dce25e237b3d5c844d306b3bb83be1))
* **export:** add declarative audit column registry with 6 sections ([a1cd1a8](https://github.com/Kohulan/ChemAudit/commit/a1cd1a8916da0b1539e5c896c8e95b79ff375f3b))
* **export:** add Excel sheet layout and SDF/PDF audit toggle to ExportDialog ([d09cf61](https://github.com/Kohulan/ChemAudit/commit/d09cf617ad3c964f8b14473ade83ac9b14dfcf65))
* **export:** add include_audit toggle to PDF report generator ([03ccdf3](https://github.com/Kohulan/ChemAudit/commit/03ccdf30eeb85133853f3e3a02f1b1dad67eac9f))
* **export:** add include_audit toggle to SDF exporter ([ff385ec](https://github.com/Kohulan/ChemAudit/commit/ff385ecf1eb3f308b4859420f30eebed40204530))
* **export:** add sheet_layout and include_audit params, update filename convention ([e98ca3f](https://github.com/Kohulan/ChemAudit/commit/e98ca3f5b19d25272816c6a5c2568fae36cfcc16))
* **export:** rewrite CSV exporter to include full audit data ([4e8713e](https://github.com/Kohulan/ChemAudit/commit/4e8713e6d8e9fd02aed0e41a26d7cc41d5d40543))
* **export:** rewrite Excel exporter with audit columns and multi-sheet mode ([498faae](https://github.com/Kohulan/ChemAudit/commit/498faae4f9bb7f8bf031ae47bd69bb9ca66f2840))
* **export:** rewrite JSON exporter with structured audit sections ([2715c2a](https://github.com/Kohulan/ChemAudit/commit/2715c2ac23fa70e7977cf3d4c9c4655170aab066))
* **export:** update ExporterFactory to forward constructor kwargs ([49e0eae](https://github.com/Kohulan/ChemAudit/commit/49e0eaec4a78bc9d3865e71c43faa0883b392e61))
* **footer:** replace inline footer with cinematic GSAP scroll-reveal footer ([77527ba](https://github.com/Kohulan/ChemAudit/commit/77527ba49b5a94ac1bc4a02f0ee7f918e91ff575))


### Bug Fixes

* **about:** remove unused ExternalLinkButton to fix CI type check ([e4f38ab](https://github.com/Kohulan/ChemAudit/commit/e4f38abd9770d16efc37b66c1d3b2aef16036655))
* **about:** remove unused ExternalLinkButton to fix CI type check ([285e3ee](https://github.com/Kohulan/ChemAudit/commit/285e3ee833c855ea3aa1064f60856eb79b4c307a))
* **batch:** add scroll guard to prevent programmatic scroll into footer ([8b2cca2](https://github.com/Kohulan/ChemAudit/commit/8b2cca2f6085d2e2f0380cb2a3f823b3c1cc9ad7))
* **deps:** resolve all npm security vulnerabilities ([17c793c](https://github.com/Kohulan/ChemAudit/commit/17c793c17311a35425645a3fa30ebcf939d9ca3a))
* **export:** align audit extractors with batch data structure, fix image placement ([929760c](https://github.com/Kohulan/ChemAudit/commit/929760c3e9a05200af69377925bca3ff9bd53694))
* resolve all ruff lint errors (unused imports, naming, sorting) ([037f627](https://github.com/Kohulan/ChemAudit/commit/037f627e473e3c34e892410e60395fa799aedd19))
* resolve all ruff lint errors (unused imports, naming, sorting) ([315688a](https://github.com/Kohulan/ChemAudit/commit/315688a7e6643db4fc7c2a6e85f8de112778bb11))
* **security:** add ownership checks to structure-filter and dataset WebSockets ([1f66b3b](https://github.com/Kohulan/ChemAudit/commit/1f66b3b294ec374faeffa97592ec50d362b3ac8d))
* **security:** add rate limit to DELETE /me/data (5/min) ([536a74b](https://github.com/Kohulan/ChemAudit/commit/536a74be4cd90baff76d7c4db7cb89aacf0d5ead))
* **security:** add rate limits to /api-keys CRUD (5-10/min) ([1cf1d17](https://github.com/Kohulan/ChemAudit/commit/1cf1d1718fcc38901a0af879602bf9e383bee550))
* **security:** add rate limits to /bookmarks endpoints (5-30/min) ([72f8f9e](https://github.com/Kohulan/ChemAudit/commit/72f8f9e7d2ea25a5af2d67246d1972b990b1f476))
* **security:** add rate limits to /config and /standardize/options (30/min) ([57aaa5d](https://github.com/Kohulan/ChemAudit/commit/57aaa5d873f15832772b3012423013966df43282))
* **security:** add rate limits to /history endpoints (30/min) ([9a61708](https://github.com/Kohulan/ChemAudit/commit/9a61708b7d312b89abb2df9d02f92dfae069354e))
* **security:** add rate limits to permalink resolve endpoints (30/min) ([fcf2b66](https://github.com/Kohulan/ChemAudit/commit/fcf2b66aefc3be07877f3d3fb65f98687a846b8e))
* **security:** prevent exception details from leaking to API responses ([a8fba76](https://github.com/Kohulan/ChemAudit/commit/a8fba76f51d5f81b767ea2f711a1335230d4f0a5))
* **security:** reduce file upload rate limit from 10/min to 3/min ([4b5f798](https://github.com/Kohulan/ChemAudit/commit/4b5f798a85aa92d579651902c577403e05b25fcc))
* **security:** stop reflecting user input in profiler error responses ([bf13574](https://github.com/Kohulan/ChemAudit/commit/bf13574d6379f2be7ef411391e8419a630a3fd95))
* **security:** stop reflecting user input in profiler error responses ([d6523ac](https://github.com/Kohulan/ChemAudit/commit/d6523ac22e99dd1a3cbd4fb55a9208bafa2bc869))
* **security:** store session ownership for structure-filter and dataset jobs ([245ca4e](https://github.com/Kohulan/ChemAudit/commit/245ca4efd878b1951c137e2bab3dae6c376e8656))
* **security:** wrap profiler service calls to prevent stack trace exposure ([8fdf39c](https://github.com/Kohulan/ChemAudit/commit/8fdf39c803793ad65f3c938b5648af49b8f3de95))
* update tests for sanitized error messages, handle missing medche… ([8361e83](https://github.com/Kohulan/ChemAudit/commit/8361e83a8a50d42c295544b35334300dc16a4757))
* update tests for sanitized error messages, handle missing medchem gracefully ([74aa273](https://github.com/Kohulan/ChemAudit/commit/74aa273ff6132c1b77faae14f1bc230ccee72d75))


### Documentation

* add SECURITY.md with vulnerability reporting policy ([abc564f](https://github.com/Kohulan/ChemAudit/commit/abc564fa8c766b39ee292feb6372c24d757a599b))
* update documentation for all new features since initial release ([be70926](https://github.com/Kohulan/ChemAudit/commit/be70926a396c0d7fd88ae592ac8133eef2dcffaf))

## [1.3.0](https://github.com/Kohulan/ChemAudit/compare/v1.2.1...v1.3.0) (2026-03-05)


### Features

* add /integrations/compare endpoint ([01649c6](https://github.com/Kohulan/ChemAudit/commit/01649c6d0abb26ae140dbd9d299420b37ce80ad5))
* add /resolve endpoint for universal identifier resolution ([d9aae97](https://github.com/Kohulan/ChemAudit/commit/d9aae97ff689e8d408eb73d4cc156560d895757c))
* add cross-database comparator service ([0a2150e](https://github.com/Kohulan/ChemAudit/commit/0a2150e62741e8aabfce78e1a35f6921a7ede205))
* add frontend types and API client methods for resolver and comparator ([d3769c2](https://github.com/Kohulan/ChemAudit/commit/d3769c2ee109936fc2ddbc29511a1ac91e37ad32))
* add identifier type detection for universal resolver ([3eea6df](https://github.com/Kohulan/ChemAudit/commit/3eea6df216d7948fd68ec3687f38dd9ed78f5f43))
* add IdentifierResolverCard and DatabaseComparisonPanel components ([0a0b4c6](https://github.com/Kohulan/ChemAudit/commit/0a0b4c65f52a6b2e95ffd03121d55caf21a2668d))
* add UniChem client for cross-database ID mapping ([30e1541](https://github.com/Kohulan/ChemAudit/commit/30e154185f24fcb217b8976e07cb4837a09fbd8f))
* add universal identifier resolver orchestrator ([b134055](https://github.com/Kohulan/ChemAudit/commit/b1340553d642bc5e94170fca6d7c96c07cc62f68))
* add Wikidata and ChEBI clients for identifier resolution ([fa3e6fc](https://github.com/Kohulan/ChemAudit/commit/fa3e6fc46b161c7842bf001bd3c532d7db645b6c))
* fix stereochemistry in cross-db comparison, add Wikidata support, update docs ([1bceb62](https://github.com/Kohulan/ChemAudit/commit/1bceb621c23b604f8c2a3296fc6065c212bac8ab))
* include dois ([f71347f](https://github.com/Kohulan/ChemAudit/commit/f71347f485c1e16331c706a00c35edc17c77edf2))
* include dois ([f1f3a23](https://github.com/Kohulan/ChemAudit/commit/f1f3a238ae48063fe3cf31042521109bc195d2c2))
* integrate comparison panel into DatabaseLookup and resolver into SingleValidation ([69d7939](https://github.com/Kohulan/ChemAudit/commit/69d7939d68c272831aac0803687f92190cd79bae))
* support multi-identifier input in single validation UI ([0f95466](https://github.com/Kohulan/ChemAudit/commit/0f95466071515ee16622eaae85bdf3d44b648bb3))


### Bug Fixes

* prevent cookie injection by never echoing user-supplied session IDs ([a0eed30](https://github.com/Kohulan/ChemAudit/commit/a0eed306bb53e6cf608d314bee7afa3423392e89))
* prevent cookie injection by never echoing user-supplied session IDs ([1b0d07d](https://github.com/Kohulan/ChemAudit/commit/1b0d07d4c4ded874e7568f23447d5cad70625933))
* resolve CodeQL security alerts for PR [#43](https://github.com/Kohulan/ChemAudit/issues/43) ([b55abbd](https://github.com/Kohulan/ChemAudit/commit/b55abbd991699f75d67d287b55887d18425f67f3))
* security fixes ([f687aac](https://github.com/Kohulan/ChemAudit/commit/f687aac7fc1796b6a9016135af8b004b7ea5c428))
* security fixes ([2e90988](https://github.com/Kohulan/ChemAudit/commit/2e90988f6c16b053a13403ab3213aee62e9e0c02))
* use PBKDF2 for API key lookup hash (CodeQL-approved KDF) ([57fc32e](https://github.com/Kohulan/ChemAudit/commit/57fc32edabbd28b8333c2a149107f4720d7a2779))

## [1.2.1](https://github.com/Kohulan/ChemAudit/compare/v1.2.0...v1.2.1) (2026-02-27)


### Bug Fixes

* allow Swagger UI CDN assets in Content-Security-Policy ([aa0cbe7](https://github.com/Kohulan/ChemAudit/commit/aa0cbe716629dac601a101b64cf1c6060fcbef79))
* allow Swagger UI CDN assets in Content-Security-Policy ([a4d69ec](https://github.com/Kohulan/ChemAudit/commit/a4d69ec51cb97acca710387c72cc58feb50d692c))
* enable API docs (Swagger, ReDoc) in production ([771524e](https://github.com/Kohulan/ChemAudit/commit/771524ed8543e60e8ae715b17b4b9a43cb9a0d9a))
* enable API docs (Swagger, ReDoc) in production ([cbba13d](https://github.com/Kohulan/ChemAudit/commit/cbba13da3dd203a6fde8f2e55b34bd9bffc0fbb3))

## [1.2.0](https://github.com/Kohulan/ChemAudit/compare/v1.1.3...v1.2.0) (2026-02-27)


### Features

* add /profiles route with PresetPicker and ProfileBuilder, nav link in header ([aae56ad](https://github.com/Kohulan/ChemAudit/commit/aae56ad7b8730d4d00f2fd91bfc418c6efa8fc32))
* add /report/:shortId permalink route for shared batch reports ([17465c2](https://github.com/Kohulan/ChemAudit/commit/17465c22f4ebb80ef36ec0433a7284c4df9b0d49))
* add /report/:shortId permalink route for shared batch reports ([2d3cc95](https://github.com/Kohulan/ChemAudit/commit/2d3cc95eb2b9b75897d1f45f2afe07c3f466e313))
* add /report/:shortId permalink route for shared batch reports ([f95b250](https://github.com/Kohulan/ChemAudit/commit/f95b250f0052024a86f966c8d88490b8fdf1e1ea))
* add Bioavailability & Permeation references to About page ([e27d29a](https://github.com/Kohulan/ChemAudit/commit/e27d29a62457e09d4ee348f40e74e41f9279f3bd))
* add bookmark result snapshots with IndexedDB storage ([86ae31a](https://github.com/Kohulan/ChemAudit/commit/86ae31abf8745ef3271a900ce8850e904c3790a3))
* add conditional Profile Score column to batch results table ([3068e1f](https://github.com/Kohulan/ChemAudit/commit/3068e1fc860354b8e96bf4db666c33e7c46c650b))
* add desirability-based profile scoring module ([5e8de57](https://github.com/Kohulan/ChemAudit/commit/5e8de5751a672f44fb3ba8d36326de1334f8c839))
* add informational tooltips to Bioavailability Radar and BOILED-Egg sections ([c396f4b](https://github.com/Kohulan/ChemAudit/commit/c396f4bc08e4babb43969487e1bfe2cd1916806d))
* add pattern category classifier for structural alerts ([f2fd38f](https://github.com/Kohulan/ChemAudit/commit/f2fd38fdeff6bb07aa8e30015d46fffed6a39a7c))
* add profile score distribution chart to batch analytics ([eddb58a](https://github.com/Kohulan/ChemAudit/commit/eddb58ae4c5624333532808b52766732ac88da18))
* add profile_id to batch upload API and types ([6784831](https://github.com/Kohulan/ChemAudit/commit/6784831146c135a17e69630c24c7c36c7f9c6af5))
* add profile_score sort extractor and profile stats to aggregator ([220d290](https://github.com/Kohulan/ChemAudit/commit/220d2902be4b1489dba5648e48db87639430a68a))
* add rotatable_bonds and aromatic_rings to batch results ([f2089da](https://github.com/Kohulan/ChemAudit/commit/f2089daefc229b6b19897dcf6ef05fc85d8043f4))
* add scoring visual polish with count-up animations and chart variants ([50af474](https://github.com/Kohulan/ChemAudit/commit/50af4741d57b4ed7fcc0927548ba5a77f3ee5e06))
* add session auto-purge, privacy page overhaul, and integration tests ([d311925](https://github.com/Kohulan/ChemAudit/commit/d3119254af24bf1d1d9ba4047991be449cf53266))
* add session-based data isolation for bookmarks and history ([9c2e7cb](https://github.com/Kohulan/ChemAudit/commit/9c2e7cbc9652b41d419b05d6ec5f0eef06a88467))
* add tooltips and result messages to all 27 validation checks ([6572263](https://github.com/Kohulan/ChemAudit/commit/6572263a4ab7658b6cd8fb6e83e4ce4dc85fed66))
* add validation cache, canonical SMILES source, and batch-to-single navigation ([ee2f934](https://github.com/Kohulan/ChemAudit/commit/ee2f9347639884252f3fa55e0cf5142dc9a1e959))
* batch UX overhaul — reorder sections, fix compare bug, add quick-nav and outlier popover ([a1a7a8b](https://github.com/Kohulan/ChemAudit/commit/a1a7a8be39298395a6841fdf46654c6c0ee93e34))
* create ProfileSidebar component for batch upload flow ([66dffee](https://github.com/Kohulan/ChemAudit/commit/66dffee79118e88bf9cda4d49e814f118a27e8c8))
* enhance batch UX with comparison panel, export dialog, and subset actions ([4b070dd](https://github.com/Kohulan/ChemAudit/commit/4b070ddd9118d4bde74dd83c5b78cee95e57d05c))
* enrich AVAILABLE_CATALOGS with citations, DOIs, and scope descriptions ([88c2c85](https://github.com/Kohulan/ChemAudit/commit/88c2c85437ec6d99635d9a8d94caecbdc0f42010))
* enrich structural alerts with catalog context, citations, and categories ([b3eba8e](https://github.com/Kohulan/ChemAudit/commit/b3eba8ef669e1df96cafc09bf41fce9cfc174c79))
* expand BOILED-Egg tooltip with model details, boundaries, and cross-validation stats ([9eb7ea3](https://github.com/Kohulan/ChemAudit/commit/9eb7ea31d40401ee690097e3d7fd6ee840adfd05))
* extract RDKit entry metadata and classify alert patterns ([2821336](https://github.com/Kohulan/ChemAudit/commit/28213369c90a85574e91e2c4ef9c1f729e2c10cc))
* implement fingerprint and deduplication exporters ([b7eec50](https://github.com/Kohulan/ChemAudit/commit/b7eec5062dc304ba109847e756daa78c24c4292c))
* implement IUPAC converter, enhanced PDF sections, and tests ([cbfd4f0](https://github.com/Kohulan/ChemAudit/commit/cbfd4f056c984b05dd9905c66c6f7f5ec128b1ad))
* implement notifications, permalinks, and audit trail ([77cd4ec](https://github.com/Kohulan/ChemAudit/commit/77cd4ecfa651fe77a250543b8507a5b8dbe323bc))
* implement profiles CRUD with 8 presets, bookmarks, and batch subset actions ([d17b4be](https://github.com/Kohulan/ChemAudit/commit/d17b4becdf57479b2af45946bc7e4e27759cf0d4))
* implement scaffold and property matrix exporters with tests ([06fcc1c](https://github.com/Kohulan/ChemAudit/commit/06fcc1cea826ff3c1947b2384639a5744c964e2b))
* improve batch analytics, chemical space scatter, and profile sidebar ([e276f43](https://github.com/Kohulan/ChemAudit/commit/e276f43b9c903e3631b347907ffc509470bbacc6))
* improve export pipeline, IUPAC conversion, and deployment config ([da9ca07](https://github.com/Kohulan/ChemAudit/commit/da9ca07d93b62d9caccb871f8c73df7c24db42b5))
* integrate Phase 6 features into frontend — profiles, bookmarks, history, IUPAC, exports ([136cd38](https://github.com/Kohulan/ChemAudit/commit/136cd38057ddc6639e795600aeb726ebe1cf8d62))
* integrate SubsetActionPanel and Share permalink button in BatchValidation ([c588ea6](https://github.com/Kohulan/ChemAudit/commit/c588ea6a6793519aba55aecf997e18f6151699dd))
* move PDF checkboxes inline below PDF card and update input placeholder ([acc9406](https://github.com/Kohulan/ChemAudit/commit/acc94067b926ea6ad75e8c0e6f16c0295bd884de))
* provision OPSIN runtime dependencies for IUPAC conversion ([89df01a](https://github.com/Kohulan/ChemAudit/commit/89df01ad39798b8ce799e6a5763893408b695f7d))
* redesign header as floating glass pill with code simplification ([4bb9216](https://github.com/Kohulan/ChemAudit/commit/4bb92169f9650866ce467df08854541a802e8e8b))
* redesign ML readiness scorer with 4-dimension scientific assessment ([17e9289](https://github.com/Kohulan/ChemAudit/commit/17e92899678857be4235f3e36637841edc5b91eb))
* restyle bookmark button with ClayButton label and clear results on tab switch ([6e6dbda](https://github.com/Kohulan/ChemAudit/commit/6e6dbda410bb635356d05f8b3145f4e9f59ca51d))
* set up SQLAlchemy ORM foundation with 4 models and Alembic ([a0b54eb](https://github.com/Kohulan/ChemAudit/commit/a0b54eb3c2d47bfeb88830f0959692a246d9271e))
* wire batch audit trail and webhook dispatch on completion ([8e74c80](https://github.com/Kohulan/ChemAudit/commit/8e74c808b278311045a7b66c40b2e48a244e6aac))
* wire email dispatch in batch aggregation tasks and add NOTIFICATION_EMAIL config ([b85ab2e](https://github.com/Kohulan/ChemAudit/commit/b85ab2eec16f0ca26bdc58100111488adfa38d5a))
* wire profile scoring into batch upload and Celery processing ([732064e](https://github.com/Kohulan/ChemAudit/commit/732064e6b24fdda81eb0f54965f72ab5fbfaebd2))
* wire ProfileSidebar into batch validation upload flow ([a69c9aa](https://github.com/Kohulan/ChemAudit/commit/a69c9aae9ea282e97c8eaff90bae391ca2b5ae7f))


### Bug Fixes

* add SETGID/SETUID capabilities to Redis container for user switc… ([84f7b45](https://github.com/Kohulan/ChemAudit/commit/84f7b456cb738291f06f02f664f8606f6a7f8eae))
* add SETGID/SETUID capabilities to Redis container for user switching ([d764874](https://github.com/Kohulan/ChemAudit/commit/d764874b8922568756e7b2880c2890d36e2db8d0))
* comprehensive security hardening across backend, frontend, and infrastructure ([8196f65](https://github.com/Kohulan/ChemAudit/commit/8196f65db20889ec020bfea0f0e83613f7168aea))
* move Actions button next to Compare and wire onNewJob to navigate to new batch job ([c1e5b35](https://github.com/Kohulan/ChemAudit/commit/c1e5b359b2108de859828e46753d0754c9a24587))
* resolve CI lint errors and test failures ([967268e](https://github.com/Kohulan/ChemAudit/commit/967268ea99b3674be84d7ffe568a804137646753))
* resolve CI lint errors and test failures ([2a40acb](https://github.com/Kohulan/ChemAudit/commit/2a40acb76237bd90dead42c4e5525ab8ae5a1053))
* resolve lint errors and clean up unused variables across backend ([cd6c1e7](https://github.com/Kohulan/ChemAudit/commit/cd6c1e79bace5e3dabf1f8bb588fb3338d58d3e3))
* resolve merge conflict in standardization types ([7c1204a](https://github.com/Kohulan/ChemAudit/commit/7c1204aa5aca22c97e7d86f425bb57b2b21e2f73))
* resolve permalink results not loading and improve scoring profile visibility ([ee57ed8](https://github.com/Kohulan/ChemAudit/commit/ee57ed882448fd60d62296bbd7362bdab088ab30))
* use set_config() instead of SET LOCAL for RLS session variables ([9daae3d](https://github.com/Kohulan/ChemAudit/commit/9daae3df364d05a7238d7d36cb0a6f2de95bd030))
* wire Alembic migrations into Docker startup and fix RLS session vars ([1f5b48b](https://github.com/Kohulan/ChemAudit/commit/1f5b48bd09f2212faf44b38d747c1df3fa4ffb90))


### Documentation

* add Bioavailability & Permeation docs page, update overview and sidebar ([6393191](https://github.com/Kohulan/ChemAudit/commit/63931916c96f8c0a9d6fe3f96672896a5398f001))
* add design and implementation plans for features_v1 work ([f424351](https://github.com/Kohulan/ChemAudit/commit/f424351b2b99ab1ac9d59129f3de52c7ebe239c2))
* add scoring methodology reference and update user guide ([ec83ca5](https://github.com/Kohulan/ChemAudit/commit/ec83ca5b0401ec51ea4d8f30d3eacb47f8cd33e7))
* complete 06-06 backend integration gap closure plan ([b412f73](https://github.com/Kohulan/ChemAudit/commit/b412f7393a440ac60c3c3b9cd75b0ad6e0a91fc4))
* complete 06-09 email dispatch wiring plan — WORK-13 done ([bebe8f6](https://github.com/Kohulan/ChemAudit/commit/bebe8f673c03a37759e83a726d67b9dc83eb46da))
* complete advanced export formats plan ([2a00492](https://github.com/Kohulan/ChemAudit/commit/2a00492858bcd82ebed6ac19f0eb778c6827bdd1))
* overhaul Docusaurus site with 5 new pages and 15 page updates ([049fa5d](https://github.com/Kohulan/ChemAudit/commit/049fa5d9c765d16a09aeb3bbb3e3a2b463247553))
* update About page with complete safety filter citations and accurate stats ([b3d36d1](https://github.com/Kohulan/ChemAudit/commit/b3d36d17a6250bce40f7ec46637a55c7a7d5f6ed))
* update privacy policy with IP logging, API key persistence, and localStorage disclosures ([9391dc3](https://github.com/Kohulan/ChemAudit/commit/9391dc3db90c2e3c60e2a9a557afacba59b84484))
* update ROADMAP and STATE for 06-08 completion — phase 6 fully complete ([0b10545](https://github.com/Kohulan/ChemAudit/commit/0b10545183a0486d18e37230d64f491afb334ab0))

## [1.1.3](https://github.com/Kohulan/ChemAudit/compare/v1.1.2...v1.1.3) (2026-02-24)


### Bug Fixes

* add missing MoleculeInfo fields to ValidationResults test mocks ([2c7cbd2](https://github.com/Kohulan/ChemAudit/commit/2c7cbd2b8894ce89d3aebbd326b3025c4a86ed3e))
* add missing MoleculeInfo fields to ValidationResults test mocks ([26b4a78](https://github.com/Kohulan/ChemAudit/commit/26b4a7848043c0f997df2cf57074320be8258b79))
* replace regex SMILES parsing with RDKit for molecule info ([ca452cd](https://github.com/Kohulan/ChemAudit/commit/ca452cdf582a9cb93256f1927b2d57026e78648d))
* replace regex SMILES parsing with RDKit for molecule info ([382c18b](https://github.com/Kohulan/ChemAudit/commit/382c18bcae8a4d82b2d895e7b8ec1ec7f960caa4))
* replace regex SMILES parsing with RDKit for molecule info ([4b7d483](https://github.com/Kohulan/ChemAudit/commit/4b7d4835239d00d87db3ba3a275fd06cf0547673))

## [1.1.2](https://github.com/Kohulan/ChemAudit/compare/v1.1.1...v1.1.2) (2026-02-12)


### Bug Fixes

* add WebSocket proxy location block for batch progress ([981cb1c](https://github.com/Kohulan/ChemAudit/commit/981cb1c1d6d3b505c8e9ee71b147031567654751))
* add WebSocket proxy location block for batch progress ([908c424](https://github.com/Kohulan/ChemAudit/commit/908c4245c4442833daab291721da385621ae9140))

## [1.1.1](https://github.com/Kohulan/ChemAudit/compare/v1.1.0...v1.1.1) (2026-02-06)


### Bug Fixes

* batch validation issue nginx ([99ff3e1](https://github.com/Kohulan/ChemAudit/commit/99ff3e14481c28fa2e818f77254480425cef909c))
* batch validation issue nginx ([b7b9f85](https://github.com/Kohulan/ChemAudit/commit/b7b9f85156103b0dd9e067c4664cfe2e0fbd7ca7))

## [1.1.0](https://github.com/Kohulan/ChemAudit/compare/v1.0.0...v1.1.0) (2026-02-06)


### Features

* include detailed references and highlighting download option ([4bb21ea](https://github.com/Kohulan/ChemAudit/commit/4bb21eaed3349ade759acd58a560367594b84f05))


### Bug Fixes

* aethetic fixes and depictions ([1fec46d](https://github.com/Kohulan/ChemAudit/commit/1fec46d550226f8e4150e862a25a510619c6aff4))
* aethetic fixes and depictions ([bea6200](https://github.com/Kohulan/ChemAudit/commit/bea6200ee9e145703fd976fbf869a109776382f3))
* docker frontend build with [roxy ([bf3e075](https://github.com/Kohulan/ChemAudit/commit/bf3e07521b5c5e9b77337fe6f3e3b62d5eeae512))
* docker frontend build with proxy ([c460c06](https://github.com/Kohulan/ChemAudit/commit/c460c06a43bf27aaba191f5577826b9490300503))
* frontend deployment issues ([4910f3e](https://github.com/Kohulan/ChemAudit/commit/4910f3e8996bcabd04eaf5f79e6408b8e280f235))
* frontend deployment issues ([a895034](https://github.com/Kohulan/ChemAudit/commit/a895034add369ef9b8da936a03579929dcd6f918))
* information on drawing molecule ([83d12aa](https://github.com/Kohulan/ChemAudit/commit/83d12aadd79485e9a6dddaf7a7e1dd3de82d91b6))

## 0.1.0 (2026-02-05)


### Features

* add batch processing frontend with real-time progress ([cecc359](https://github.com/Kohulan/ChemAudit/commit/cecc359911b4d131aa2c8b22d8664d0a4bb023b1))
* add Celery batch processing infrastructure ([84d6a6d](https://github.com/Kohulan/ChemAudit/commit/84d6a6d0c153bd074afb74cb0519eb8882ac6833))
* add CI pipeline and test infrastructure ([5b1374f](https://github.com/Kohulan/ChemAudit/commit/5b1374f3e61923e0ee3f4427f5d032688a3e5d78))
* add client examples and tests ([1c0137b](https://github.com/Kohulan/ChemAudit/commit/1c0137b19eace2c72a9b4e79de0b51c4874a97d3))
* add common schemas and health endpoint ([2d28580](https://github.com/Kohulan/ChemAudit/commit/2d28580e2dfb09d774d778dd716597acfa663563))
* add Ctrl+Enter keyboard shortcut with react-hotkeys-hook ([5a29625](https://github.com/Kohulan/ChemAudit/commit/5a29625c8cb7b209a67934a928dfea2d7d4a765b))
* add DECIMER and COCONUT integration clients ([004d06a](https://github.com/Kohulan/ChemAudit/commit/004d06ada8dfea1a45f39c20545d929b1c4d45a8))
* add environment configuration files ([8d54150](https://github.com/Kohulan/ChemAudit/commit/8d541507bb533b3c576b24976f00497cf502b889))
* add export API endpoint and tests ([3a93087](https://github.com/Kohulan/ChemAudit/commit/3a930873b9bbab6504c62ac79364e964f7de1d6c))
* add export dialog and PDF API integration ([6527759](https://github.com/Kohulan/ChemAudit/commit/65277597bb7ad4f4bd8a371432db7b50012aaf63))
* add indices filtering to export endpoint (GET and POST) ([633b2a4](https://github.com/Kohulan/ChemAudit/commit/633b2a4323692a0f2d99301a791585d4aadda5a2))
* add integration API endpoints with rate limiting ([65543a5](https://github.com/Kohulan/ChemAudit/commit/65543a5d4c0ceafff01e61b1841d957b132ab9a7))
* add ML-readiness calculation explanation ([e32747b](https://github.com/Kohulan/ChemAudit/commit/e32747b1774598861ac645024ec42316a77fde52))
* add ML-readiness scoring, NP-likeness, and scaffold extraction ([0ea0227](https://github.com/Kohulan/ChemAudit/commit/0ea02278024810bb26c6461ecd0dc5eb41bc6a8a))
* add Nginx reverse proxy configuration ([c1e324a](https://github.com/Kohulan/ChemAudit/commit/c1e324a733574abdb33cbd3a7481b5eb5a3d744f))
* add PDF report generator with WeasyPrint ([7d9bcc6](https://github.com/Kohulan/ChemAudit/commit/7d9bcc6e6afae202f8997b80265ee7f2d3eab020))
* add production docker-compose and deployment guide ([8229e9c](https://github.com/Kohulan/ChemAudit/commit/8229e9c49587943693fffdc42ad6c2300f5b7bbf))
* add production Dockerfiles with multi-stage builds ([fe16dfc](https://github.com/Kohulan/ChemAudit/commit/fe16dfc81da441f55af29baf0acc8f7f6def1578))
* add profiling infrastructure and baseline benchmarks ([4526c25](https://github.com/Kohulan/ChemAudit/commit/4526c2578171f9018ac81165f1e28a430274daea))
* add Prometheus metrics to FastAPI backend ([82503d3](https://github.com/Kohulan/ChemAudit/commit/82503d34301785f4c94bfbc982a2230f82fb4346))
* add PubChem and ChEMBL integrations with tests ([f6f0d83](https://github.com/Kohulan/ChemAudit/commit/f6f0d83c0e161ec63737dd29cd282f6eade7a0d3))
* add RecentMolecules dropdown component ([c49637e](https://github.com/Kohulan/ChemAudit/commit/c49637ef37d08acdc7c4bcf9718839d0236f8ed7))
* add Redis validation cache with InChIKey keys ([d9bd86a](https://github.com/Kohulan/ChemAudit/commit/d9bd86af7ed246b4c5be032460e4cffeec176a02))
* add scoring display frontend components ([859babb](https://github.com/Kohulan/ChemAudit/commit/859babb2045236a8d5d9afa23c9967e50df418f1))
* add selection badge, Export Selected button, and indices to export ([de3f530](https://github.com/Kohulan/ChemAudit/commit/de3f5303ac8ce00ebc62d108155aae500c247521))
* add selection state and checkbox column to batch results table ([54557ad](https://github.com/Kohulan/ChemAudit/commit/54557add4c1be98c0257c17439ab16fa42b487bf))
* add standardization frontend components ([931dc6c](https://github.com/Kohulan/ChemAudit/commit/931dc6c0ab8867b43ea7a633313b05dadcb9348a))
* add structural alert frontend components ([dbb7802](https://github.com/Kohulan/ChemAudit/commit/dbb780267da7b8c28abfa85f955acd8cf85e00f0))
* add structural alert screening backend ([48624f8](https://github.com/Kohulan/ChemAudit/commit/48624f8a4d2a16ea2b848b88e496a7511360b0e4))
* add TypeScript types for environment variables ([dc2cdca](https://github.com/Kohulan/ChemAudit/commit/dc2cdca1ee85eeb2a2a2c80de850687d68f71e3c))
* add URL sharing with useSearchParams ([e3c5f5a](https://github.com/Kohulan/ChemAudit/commit/e3c5f5a9bc6099982d7f3713926aa7dd6a87ef61))
* add useLocalStorage and useRecentMolecules hooks ([3009cdb](https://github.com/Kohulan/ChemAudit/commit/3009cdb64a088c09f2e0d964f18a15a5112c58ed))
* add WebSocket progress and batch API endpoints ([4a3c63f](https://github.com/Kohulan/ChemAudit/commit/4a3c63f5914d889725569eb4ef190919c61cf0a3))
* Batch processing supports tsv and txt files ([83f84e1](https://github.com/Kohulan/ChemAudit/commit/83f84e1bd50e022519eb782f0d216a7561766869))
* brand ChemAudit with crimson "Chem" across docs site ([a7d068f](https://github.com/Kohulan/ChemAudit/commit/a7d068f580bf11900defff1557e8541075f7e146))
* brand ChemAudit with crimson "Chem" across docs site ([308122b](https://github.com/Kohulan/ChemAudit/commit/308122b7018af597b5bffc90d48b7f32335092cb))
* configure Prometheus and Grafana services ([2eabfff](https://github.com/Kohulan/ChemAudit/commit/2eabfff8c025259ac57731517cb31dc26860a685))
* create API client and TypeScript types ([2491767](https://github.com/Kohulan/ChemAudit/commit/2491767b32a98c0463624339404768ff4ea5ecc1))
* create chemistry-inspired color palette and visual refresh ([d9c6c67](https://github.com/Kohulan/ChemAudit/commit/d9c6c67101e31af9084a46d1be5c462dc01a9bcc))
* create claymorphism landing page ([557bbba](https://github.com/Kohulan/ChemAudit/commit/557bbba4fb9203a53c895d68ee7ceae8962c2143))
* create claymorphism landing page ([f1c7eb6](https://github.com/Kohulan/ChemAudit/commit/f1c7eb657e9291ce305b0fa10c484b94925e92c6))
* create Docker Compose and backend scaffolding ([69d4d0b](https://github.com/Kohulan/ChemAudit/commit/69d4d0b46351c9ac6625b9e3a4c4ee48de77dff2))
* create enhanced score visualizations with informative tooltips ([1880802](https://github.com/Kohulan/ChemAudit/commit/1880802654bdc643589bfa49aa02b960bf2e1915))
* create ErrorFallback component and NotFound page ([148a401](https://github.com/Kohulan/ChemAudit/commit/148a4019d318ca4ba64ad94c1a1ddbc8c4972d0c))
* create export service with factory pattern ([edd0ac7](https://github.com/Kohulan/ChemAudit/commit/edd0ac7face8a8d9dc2b4fcda2ec8a30e342f19d))
* create Grafana monitoring dashboard ([1112739](https://github.com/Kohulan/ChemAudit/commit/11127393b51529a60a755c9baed7c3938a13aa96))
* create molecule parser with defensive sanitization ([4853b57](https://github.com/Kohulan/ChemAudit/commit/4853b57393868eaa85b1d43eaddac113340d65dd))
* create RDKit.js integration with WASM memory management ([c805b33](https://github.com/Kohulan/ChemAudit/commit/c805b33a6d7bdd53fd306c7068db50fdb65a0f81))
* create React + TypeScript + Vite frontend scaffolding ([9f0b422](https://github.com/Kohulan/ChemAudit/commit/9f0b4221fdd1d67f4a57d64c80c0288751fef073))
* create single validation page and complete workflow ([85c987f](https://github.com/Kohulan/ChemAudit/commit/85c987f2327d541e466f22e55d608ca230d4e8fa))
* create theme CSS system and sidebar navigation ([1adcc21](https://github.com/Kohulan/ChemAudit/commit/1adcc2131bb69ba288f98b9b8b5fbbfb9d981536))
* create theme CSS system and sidebar navigation ([b3a9c86](https://github.com/Kohulan/ChemAudit/commit/b3a9c861747230092903330a112b96fc56adad32))
* create UI components and wire up app ([d6faf67](https://github.com/Kohulan/ChemAudit/commit/d6faf6744e52f5fb1203026826257df90a1da4e8))
* create validation results UI components ([001c006](https://github.com/Kohulan/ChemAudit/commit/001c006f90800eee5696235c1998325a04aa09a3))
* Detailed implmentation of batch processing with user desired selections ([f07fdb9](https://github.com/Kohulan/ChemAudit/commit/f07fdb91c4b20f47638f62e9261e6f96417d186e))
* documents added as github pages ([415a58e](https://github.com/Kohulan/ChemAudit/commit/415a58ed777f241bc03beafd4a61e1545810e209))
* documents added as github pages ([d5c5891](https://github.com/Kohulan/ChemAudit/commit/d5c58917ee89621e35140abdabdf990bb3ec60fb))
* implement API key authentication and apply rate limits ([b2f37c8](https://github.com/Kohulan/ChemAudit/commit/b2f37c893f9d220f7f79b7b0c5d5d46f411160d8))
* implement lazy-loaded routes with code splitting ([8ecaaa7](https://github.com/Kohulan/ChemAudit/commit/8ecaaa7cd32fcf5fee7f9bd7b5480ca21c928c3b))
* implement SlowAPI rate limiter with Redis storage ([359bf46](https://github.com/Kohulan/ChemAudit/commit/359bf46c200b670492ada05a19aafa9e457226b2))
* implement stereochemistry and representation checks ([196dfad](https://github.com/Kohulan/ChemAudit/commit/196dfada6c4a03c0fb0b3b0f7ac62a26e77d7214))
* increase batch limits to 1GB/1M molecules and extend timeouts ([ce48071](https://github.com/Kohulan/ChemAudit/commit/ce48071b603cbdc075fa13e98d6057e54bf808b2))
* integrate ErrorBoundary and 404 route in App ([bb5cd12](https://github.com/Kohulan/ChemAudit/commit/bb5cd12ea5ef01eb7a23c588f3649f024e0bb000))
* integrate recent molecules in SingleValidation page ([16f6d7e](https://github.com/Kohulan/ChemAudit/commit/16f6d7ed855f76bdbb5b883314db5eff77a5b647))
* major improvements and bug fixes ([f37fca5](https://github.com/Kohulan/ChemAudit/commit/f37fca54e82f913b91cbc8195ed0957efb3870eb))
* matomo intergration ([5e730b2](https://github.com/Kohulan/ChemAudit/commit/5e730b2c3e3624d9a2edd954d364dc42a235b5a5))
* readme updates ([ee4116b](https://github.com/Kohulan/ChemAudit/commit/ee4116b732205813d49c08167df128e1f91b1d9c))
* scaffold Docusaurus and configure for GitHub Pages ([7c7d934](https://github.com/Kohulan/ChemAudit/commit/7c7d9346b182fc6f486827308f35a4c66376c3c1))
* scaffold Docusaurus and configure for GitHub Pages ([8ea0eb7](https://github.com/Kohulan/ChemAudit/commit/8ea0eb76e084abf8d262b03904eebe2f068bdc94))
* **scoring:** enhance ML-readiness scorer with 7 fingerprints and 451 descriptors ([bc1ea14](https://github.com/Kohulan/ChemAudit/commit/bc1ea1424204826f07212c16040ab0b175488eef))
* update API service with environment configuration ([2963b32](https://github.com/Kohulan/ChemAudit/commit/2963b320eb465bb30dc178aae5c2a6912ee01a10))
* update screenshots ([e315e70](https://github.com/Kohulan/ChemAudit/commit/e315e70e739eb60d237f9e881d78604bdf437149))
* use ChemAudit logo as docs site favicon ([854fe8e](https://github.com/Kohulan/ChemAudit/commit/854fe8e887153a89759169dc9adbe021d5da6b9c))
* use ChemAudit logo as docs site favicon ([2461127](https://github.com/Kohulan/ChemAudit/commit/24611277d17f7009ffa4be780d9252982547c712))
* verify SMARTS precompilation and document final benchmarks ([1e679c9](https://github.com/Kohulan/ChemAudit/commit/1e679c9343940e47ab7a935189cd2428a3293261))
* wire up main app with routes and exception handlers ([267370c](https://github.com/Kohulan/ChemAudit/commit/267370c2744f5dc6ed7d71858bb5735895444246))


### Bug Fixes

* add GitHub logo to navbar and center-align footer ([37ee5d3](https://github.com/Kohulan/ChemAudit/commit/37ee5d372dd6a1b3e6c38ede7afdc376ce2dada2))
* Add token to release workflow for authentication ([be7730d](https://github.com/Kohulan/ChemAudit/commit/be7730d9351052719790bd2dbb258dc2608a7240))
* align footer with Docs left, Community right, copyright centered ([6794485](https://github.com/Kohulan/ChemAudit/commit/6794485bc5a7ce271aad0b5f7c7d7cdaf838c660))
* Backend linting issues ([5d1c164](https://github.com/Kohulan/ChemAudit/commit/5d1c16459c51167621f94175e7e9568b93434bd1))
* **backend:** correct Murcko scaffold generic extraction ([9a60cfd](https://github.com/Kohulan/ChemAudit/commit/9a60cfd043193f6dc6404ada08b47d189984a97a))
* cleanup frontend to adhere to linting warnings ([a777725](https://github.com/Kohulan/ChemAudit/commit/a777725c67f06ee476a1ecbc5ae8e7f6d2dd843f))
* cleanup unused frontend functions ([f80bd8a](https://github.com/Kohulan/ChemAudit/commit/f80bd8a490b8bf321ad482555abf23afad98c732))
* COCONUT tests ([6344276](https://github.com/Kohulan/ChemAudit/commit/63442769e94d865ae85af4071e9fe9f3795e762d))
* correct ErrorFallback type compatibility with react-error-boundary ([9b4f788](https://github.com/Kohulan/ChemAudit/commit/9b4f788ef7c45c6f4b32152756b801a227a2a268))
* correct RDKit descriptor count from 217 to ~208 ([5098b90](https://github.com/Kohulan/ChemAudit/commit/5098b90d411fb80369b4f08f541eb452a513776a))
* docker compose production ([c011659](https://github.com/Kohulan/ChemAudit/commit/c0116593c2f3f289518c8343c19cddf01bbb5d88))
* docker production deployment ([a26258d](https://github.com/Kohulan/ChemAudit/commit/a26258d4df53afb8230395607965a641793d4761))
* **docker:** correct API URL and healthcheck endpoint ([c06f22c](https://github.com/Kohulan/ChemAudit/commit/c06f22cbf882cff8d88575c5c0140d591f4f7471))
* docs footer with coffee icon, linked author and Steinbeck Lab ([e36499e](https://github.com/Kohulan/ChemAudit/commit/e36499e3cb0d0864638a0e776d5b25f583d102ac))
* docs footer with coffee icon, linked author and Steinbeck Lab ([52bc23f](https://github.com/Kohulan/ChemAudit/commit/52bc23f0051be5edbdcc83a64ef1a7f8132d9881))
* export tests ([fe18ead](https://github.com/Kohulan/ChemAudit/commit/fe18ead13a56f865b869807e6e538384d43db440))
* **frontend:** link API Docs to backend /docs endpoint ([814e2e4](https://github.com/Kohulan/ChemAudit/commit/814e2e45e737ddba8537af9c86ef8c72f8b2e0aa))
* **frontend:** molecule sizing and error messages ([7ddf3c2](https://github.com/Kohulan/ChemAudit/commit/7ddf3c26f8c7b4df27fcc060bb15b92ec895bff3))
* make all batch results table columns sortable ([9ff5da9](https://github.com/Kohulan/ChemAudit/commit/9ff5da9a2a1c41c5368e619588eca6c0e7eda62c))
* pagination issues ([916438c](https://github.com/Kohulan/ChemAudit/commit/916438c1a2405e3ab229d6a2ee9100bc15e71696))
* project production deployment setup ([e2b4d43](https://github.com/Kohulan/ChemAudit/commit/e2b4d4334b8669582d701028e42f4b6305de1c10))
* register stereo and representation checks in validation engine ([604c8e1](https://github.com/Kohulan/ChemAudit/commit/604c8e123fa31dcbfdbd9d9fdcdf5764b4db40a6))
* resolve CI failures - add lib files, fix ruff, update workflows ([d5f32d4](https://github.com/Kohulan/ChemAudit/commit/d5f32d45a12c74dc1b2fbb6c50948f7a801865f8))
* resolve test failures and add pytest-benchmark support ([fc2afcd](https://github.com/Kohulan/ChemAudit/commit/fc2afcdddd6831c16ac65575c5a315b69b86a7d6))
* safety tests ([4446d0d](https://github.com/Kohulan/ChemAudit/commit/4446d0d7e22975d6059b0d71b8de0a0c86379b69))
* screenshots ([30878c5](https://github.com/Kohulan/ChemAudit/commit/30878c57e543ce95be3221e25e26aa12bdbbb280))
* update docs site baseUrl and links for ChemAudit rename ([6e0354d](https://github.com/Kohulan/ChemAudit/commit/6e0354d46abc4873881d2a8f0fbf0e82a27a2895))
* update docs site footer to match main website style ([8391f07](https://github.com/Kohulan/ChemAudit/commit/8391f07ee1d71d73ba776bbf0dbdb9bc7d5fafea))
* update docs site footer to match main website style ([9797d23](https://github.com/Kohulan/ChemAudit/commit/9797d23261b5bffffa23386ee85f9ac88d025b79))
* update GitHub URLs from placeholder to Kohulan/ChemAudit ([e9581c4](https://github.com/Kohulan/ChemAudit/commit/e9581c48eba31d033ad9a9eb52e6c7ccc3744925))
* update RDKit version to match backend ([14a5c61](https://github.com/Kohulan/ChemAudit/commit/14a5c6141f40d7cd9a33ed543b690b773470e815))
* use API_DOCS_URL environment variable in Header ([3b83d01](https://github.com/Kohulan/ChemAudit/commit/3b83d0139b20f95020444bc50ce4e8df3dd79517))
* use original database logos ([deaf2b1](https://github.com/Kohulan/ChemAudit/commit/deaf2b166720b2e6fee099bc8801733941250129))


### Documentation

* add comprehensive README with quick start guide ([157b9b5](https://github.com/Kohulan/ChemAudit/commit/157b9b5f1151720ed0dced317f4405af9962bddc))
* convert API Reference, Deployment, and Troubleshooting documentation ([2a4ab25](https://github.com/Kohulan/ChemAudit/commit/2a4ab25188b0543ba709e0086a7f3c87b4d00d17))
* convert API Reference, Deployment, and Troubleshooting documentation ([2439201](https://github.com/Kohulan/ChemAudit/commit/243920146ab751352028a84d4805e8c5b6e32be8))
* convert Getting Started, Introduction, and User Guide documentation to Docusaurus ([da58ba7](https://github.com/Kohulan/ChemAudit/commit/da58ba77e5f8525251bd70149a796ee9fece3a49))
* convert Getting Started, Introduction, and User Guide documentation to Docusaurus ([74c1ad8](https://github.com/Kohulan/ChemAudit/commit/74c1ad8f9f322df03744993e1326be64f02a7091))
* create comprehensive troubleshooting FAQ ([3865f96](https://github.com/Kohulan/ChemAudit/commit/3865f96887da60199feca2b368d29dcd8d8876be))
* create comprehensive user guide documentation ([01d0ba5](https://github.com/Kohulan/ChemAudit/commit/01d0ba5536a2dfe1e8b56d2e6d4fb1b022578332))
* create deployment documentation ([88627cd](https://github.com/Kohulan/ChemAudit/commit/88627cdff5f2bb0952594db365f2276b21c12b98))
* initialize project ([8a3a026](https://github.com/Kohulan/ChemAudit/commit/8a3a026593a248ecc8b85bb8f37d255200e74457))
* map existing codebase ([518120d](https://github.com/Kohulan/ChemAudit/commit/518120de4bbf8d24724c45bff1c0dd9c3162c1ff))
* update README footer attribution ([a1c6b1e](https://github.com/Kohulan/ChemAudit/commit/a1c6b1ec07664c55e35a889dfc90c54c33d117b8))
