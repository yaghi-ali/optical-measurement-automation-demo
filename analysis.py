"""Batch Gaussian beam-profile fitting, CSV ingestion and residual reporting."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit


def gaussian(x, amplitude, center, sigma, offset):
    return offset+amplitude*np.exp(-0.5*((x-center)/sigma)**2)


def fit_profile(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.ndim != 1 or x.shape != y.shape or len(x) < 10:
        raise ValueError('At least ten paired 1D samples required')
    if not np.isfinite(x).all() or not np.isfinite(y).all() or np.any(np.diff(x) <= 0):
        raise ValueError('Finite values and strictly increasing positions required')
    if np.ptp(y) <= 1e-10:
        raise ValueError('Flat signal has no identifiable beam')
    span = np.ptp(x)
    p0 = [np.ptp(y), x[np.argmax(y)], span/8, np.min(y)]
    params, cov = curve_fit(gaussian, x, y, p0=p0,
                           bounds=([0, x[0], span/10000, -np.inf],
                                   [np.inf, x[-1], span, np.inf]), maxfev=10000)
    fit = gaussian(x, *params)
    fwhm = 2*np.sqrt(2*np.log(2))*params[2]
    return dict(amplitude=float(params[0]), center_mm=float(params[1]),
                sigma_mm=float(params[2]), offset=float(params[3]), fwhm_mm=float(fwhm),
                residual_rmse=float(np.sqrt(np.mean((y-fit)**2))),
                fwhm_fit_standard_error_mm=float(2*np.sqrt(2*np.log(2))*np.sqrt(cov[2, 2]))), fit


def generate_csv(path, seed=42):
    rng = np.random.default_rng(seed)
    with Path(path).open('w', newline='') as f:
        writer = csv.writer(f); writer.writerow(['run_id', 'position_mm', 'intensity_au'])
        for i, sigma in enumerate([0.40, 0.55, 0.70]):
            x = np.linspace(-3, 3, 201)
            y = gaussian(x, 1, 0.05*i, sigma, 0.08)+rng.normal(0, 0.01, x.shape)
            writer.writerows((f'run_{i+1}', float(a), float(b)) for a, b in zip(x, y))


def read_csv(path):
    runs = {}
    with Path(path).open(newline='') as f:
        reader = csv.DictReader(f)
        required = {'run_id', 'position_mm', 'intensity_au'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Required CSV columns: '+', '.join(sorted(required)))
        for row in reader:
            runs.setdefault(row['run_id'], []).append((float(row['position_mm']), float(row['intensity_au'])))
    if not runs:
        raise ValueError('Empty measurement file')
    return {key: np.array(sorted(rows)) for key, rows in runs.items()}


def run(out='results', source=None, seed=42):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    if source is None:
        source = out/'synthetic_measurements.csv'; generate_csv(source, seed)
        provenance = 'synthetic Gaussian profiles with additive Gaussian noise'
    else:
        provenance = 'user-supplied CSV; provenance not independently verified'
    runs = read_csv(source)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    results = []
    for name, data in runs.items():
        metrics, fitted = fit_profile(data[:, 0], data[:, 1])
        results.append(dict(run_id=name, **metrics))
        ax[0].plot(data[:, 0], data[:, 1], '.', alpha=0.25)
        ax[0].plot(data[:, 0], fitted, label=name)
        ax[1].plot(data[:, 0], data[:, 1]-fitted, label=name)
    ax[0].set(xlabel='Position (mm)', ylabel='Intensity (a.u.)', title='Gaussian fits'); ax[0].legend()
    ax[1].set(xlabel='Position (mm)', ylabel='Residual (a.u.)', title='Model residuals')
    fig.suptitle('Optical measurement analysis | '+('synthetic demo' if provenance.startswith('synthetic') else 'CSV input'))
    fig.savefig(out/'preview.png', dpi=150); plt.close(fig)
    report = dict(data=provenance, seed=seed if provenance.startswith('synthetic') else None, runs=results)
    (out/'metrics.json').write_text(json.dumps(report, indent=2))
    with (out/'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys()); writer.writeheader(); writer.writerows(results)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='results')
    parser.add_argument('--csv', default=None)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    print(json.dumps(run(args.out, args.csv, args.seed), indent=2))
