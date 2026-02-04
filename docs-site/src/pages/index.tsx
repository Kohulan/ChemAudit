import React from 'react';
import Layout from '@theme/Layout';

export default function Home(): JSX.Element {
  return (
    <Layout
      title="ChemAudit"
      description="Chemical Structure Validation Suite">
      <main>
        <div style={{ padding: '2rem', textAlign: 'center' }}>
          <h1>ChemAudit</h1>
          <p>Chemical Structure Validation Suite</p>
          <p>Documentation coming soon.</p>
        </div>
      </main>
    </Layout>
  );
}
