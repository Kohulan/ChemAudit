# Changelog

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
