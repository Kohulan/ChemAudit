import { describe, expect, it } from 'vitest';
import { csvCell } from './csv';

describe('csvCell', () => {
  it('prefixes text a spreadsheet would run as a formula', () => {
    expect(csvCell('=1+1')).toBe("'=1+1");
    expect(csvCell('@SUM(A1)')).toBe("'@SUM(A1)");
    expect(csvCell('-A1')).toBe("'-A1");
    expect(csvCell('{=A1}')).toBe("'{=A1}");
  });

  it('leaves numbers and normal text alone', () => {
    expect(csvCell('-1.5')).toBe('-1.5');
    expect(csvCell('CCO')).toBe('CCO');
    expect(csvCell(-2)).toBe('-2');
    expect(csvCell(true)).toBe('true');
  });

  it('returns an empty string for null and undefined', () => {
    expect(csvCell(null)).toBe('');
    expect(csvCell(undefined)).toBe('');
  });

  it('quotes values with commas, quotes or newlines', () => {
    expect(csvCell('a,b')).toBe('"a,b"');
    expect(csvCell('say "hi"')).toBe('"say ""hi"""');
    expect(csvCell('=a,b')).toBe('"\'=a,b"');
    expect(csvCell('a\nb')).toBe('"a\nb"');
  });
});
