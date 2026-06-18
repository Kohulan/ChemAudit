import {themes as prismThemes} from 'prism-react-renderer';
import type {Config} from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';
import matter from 'gray-matter';
import yaml from 'js-yaml';

const config: Config = {
  title: 'ChemAudit',
  tagline: 'Chemical Structure Validation Suite',
  favicon: 'img/favicon.ico',

  // Set the production url of your site here
  url: 'https://kohulan.github.io',
  // Set the /<baseUrl>/ pathname under which your site is served
  // For GitHub pages deployment, it is often '/<projectName>/'
  baseUrl: '/ChemAudit/',

  // GitHub pages deployment config.
  organizationName: 'Kohulan',
  projectName: 'ChemAudit',
  trailingSlash: false,

  onBrokenLinks: 'throw',

  markdown: {
    hooks: {
      onBrokenMarkdownLinks: 'warn',
    },
    // Docusaurus' default front matter parser calls gray-matter, whose default
    // YAML engine uses yaml.safeLoad — removed in js-yaml 4. We pin js-yaml to
    // the patched 4.x (GHSA-h67p-54hq-rp68) everywhere, so we supply a js-yaml 4
    // engine here. Mirrors DEFAULT_PARSE_FRONT_MATTER (structuredClone + trim).
    parseFrontMatter: async ({fileContent}) => {
      const {data, content} = matter(fileContent, {
        engines: {
          yaml: {
            parse: (str: string) => yaml.load(str) as object,
            stringify: (obj: object) => yaml.dump(obj),
          },
        },
      });
      return {frontMatter: structuredClone(data), content: content.trim()};
    },
  },

  // Even if you don't use internationalization, you can use this field to set
  // useful metadata like html lang. For example, if your site is Chinese, you
  // may want to replace "en" with "zh-Hans".
  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  stylesheets: [
    {
      href: 'https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap',
      type: 'text/css',
    },
  ],

  presets: [
    [
      'classic',
      {
        docs: {
          sidebarPath: './sidebars.ts',
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      } satisfies Preset.Options,
    ],
  ],

  themeConfig: {
    // Replace with your project's social card
    image: 'img/docusaurus-social-card.jpg',
    colorMode: {
      defaultMode: 'light',
      disableSwitch: false,
      respectPrefersColorScheme: true,
    },
    docs: {
      sidebar: {
        hideable: true,
      },
    },
    navbar: {
      title: '',
      logo: {
        alt: 'ChemAudit Logo',
        src: 'img/logo.png',
      },
      items: [
        {
          type: 'docSidebar',
          sidebarId: 'docs',
          position: 'left',
          label: 'Docs',
        },
        {
          href: 'https://github.com/Kohulan/ChemAudit',
          position: 'right',
          className: 'header-github-link',
          'aria-label': 'GitHub repository',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Getting Started',
              to: '/docs/getting-started/installation',
            },
            {
              label: 'User Guide',
              to: '/docs/user-guide/single-validation',
            },
            {
              label: 'API Reference',
              to: '/docs/api/overview',
            },
          ],
        },
        {
          title: 'Community',
          items: [
            {
              label: 'GitHub',
              href: 'https://github.com/Kohulan/ChemAudit',
            },
          ],
        },
      ],
      copyright: `Made with ☕ by <a href="https://kohulanr.com" target="_blank" rel="noopener noreferrer">Kohulan.R</a> at <a href="http://cheminf.uni-jena.de/" target="_blank" rel="noopener noreferrer">Steinbeck Lab</a> · © ${new Date().getFullYear()} <a href="https://chemaudit.naturalproducts.net" target="_blank" rel="noopener noreferrer"><span style="color:#E53456;font-weight:800">Chem</span>Audit</a>`,
    },
    prism: {
      theme: prismThemes.github,
      darkTheme: prismThemes.dracula,
    },
  } satisfies Preset.ThemeConfig,
};

export default config;
