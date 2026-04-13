# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 1.4.x   | :white_check_mark: |
| 1.3.x   | :white_check_mark: |
| < 1.3   | :x:                |

Only the two most recent minor releases receive security patches.
Older versions should be upgraded.

## Reporting a Vulnerability

**Please do not open a public GitHub issue for security vulnerabilities.**

If you discover a security issue in ChemAudit, report it privately using
one of the following methods:

1. **GitHub Private Vulnerability Reporting** (preferred)
   Navigate to the [Security Advisories](https://github.com/Kohulan/ChemAudit/security/advisories/new)
   page and click **"Report a vulnerability"**.

2. **Email**
   Send details to **kohulan.rajan@uni-jena.de** with the subject line
   `[ChemAudit Security]`.

Include the following in your report:

- Description of the vulnerability
- Steps to reproduce
- Affected version(s)
- Potential impact assessment
- Any suggested fix (optional but appreciated)

## Response Timeline

| Stage               | Target     |
| ------------------- | ---------- |
| Acknowledgment      | 48 hours   |
| Initial assessment  | 5 business days |
| Patch release       | 30 days (critical: 7 days) |

We will keep you informed of progress and coordinate disclosure timing
with you before any public announcement.

## Scope

The following are considered in-scope for security reports:

- Authentication and authorization bypass (API key handling, session management)
- Injection vulnerabilities (SQL, command, SMILES-based injection)
- Cross-site scripting (XSS) or cross-site request forgery (CSRF)
- Sensitive data exposure (API keys, user data, batch results)
- Rate limiting or IP ban bypass
- Insecure deserialization
- Server-side request forgery (SSRF)
- Dependency vulnerabilities with a demonstrated exploit path

The following are **out of scope**:

- Denial of service via high-volume requests (use rate limiting configuration)
- Issues in third-party dependencies without a proof of concept against ChemAudit
- Social engineering
- Physical security

## Security Architecture

ChemAudit implements the following security controls:

- **API key authentication** with Redis-backed storage and rotation support
- **CSRF protection** on state-changing endpoints
- **Rate limiting** (slowapi) with configurable thresholds
- **IP banning** for repeated abuse
- **Input validation** via Pydantic v2 schemas on all API boundaries
- **Dependency scanning** via GitHub Dependabot

## Acknowledgments

We appreciate the security research community's efforts in helping keep
ChemAudit and its users safe. Reporters who follow responsible disclosure
will be credited in release notes (unless they prefer to remain anonymous).
