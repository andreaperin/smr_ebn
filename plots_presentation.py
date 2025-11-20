"""Visualize probabilistic sampling versus uniform coverage with separate plots."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def sample_from_distributions(num_points: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Draw x and y samples from two independent non-uniform distributions."""
    x = rng.normal(loc=0.0, scale=1.0, size=num_points)
    # Beta distribution keeps most mass near 0, so corners remain under-sampled.
    y = rng.beta(a=2.0, b=5.0, size=num_points) * 6 - 3  # shift to match x-axis range
    return x, y


def uniform_grid(num_ticks: int, x_limits: tuple[float, float], y_limits: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    """Generate meshgrid coordinates that cover the full rectangle uniformly."""
    x_vals = np.linspace(*x_limits, num_ticks)
    y_vals = np.linspace(*y_limits, num_ticks)
    grid_x, grid_y = np.meshgrid(x_vals, y_vals)
    return grid_x.ravel(), grid_y.ravel()


def response_surface(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Evaluate a smooth function that depends on both inputs."""
    return np.sin(x) * np.exp(-0.3 * y**2)


def main() -> None:
    rng = np.random.default_rng(seed=3)
    num_samples = 600

    x_samples, y_samples = sample_from_distributions(num_samples, rng)
    x_limits = (-3, 3)
    y_limits = (-3, 3)

    x_uniform, y_uniform = uniform_grid(num_ticks=25, x_limits=x_limits, y_limits=y_limits)

    z_samples = response_surface(x_samples, y_samples)
    z_uniform = response_surface(x_uniform, y_uniform)

    fig_x = plt.figure(figsize=(6, 4))
    ax_x_hist = fig_x.add_subplot(111)
    ax_x_hist.hist(x_samples, bins=30, density=True, alpha=0.7, color="tab:blue")
    ax_x_hist.set_xlabel("x")
    ax_x_hist.set_ylabel("density")
    fig_x.tight_layout()

    fig_y = plt.figure(figsize=(6, 4))
    ax_y_hist = fig_y.add_subplot(111)
    ax_y_hist.hist(y_samples, bins=30, density=True, alpha=0.7, color="tab:green")
    ax_y_hist.set_xlabel("y")
    ax_y_hist.set_ylabel("density")
    fig_y.tight_layout()

    fig_sample = plt.figure(figsize=(7, 5))
    ax_sample = fig_sample.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("viridis")
    scatter_sample = ax_sample.scatter(
        x_samples,
        y_samples,
        z_samples,
        c=z_samples,
        cmap=cmap,
        alpha=0.7,
        s=20,
        edgecolor="white",
        linewidth=0.2,
    )
    ax_sample.set_xlabel("x")
    ax_sample.set_ylabel("y")
    ax_sample.set_zlabel("f(x, y)")
    ax_sample.set_xlim(*x_limits)
    ax_sample.set_ylim(*y_limits)
    fig_sample.colorbar(scatter_sample, ax=ax_sample, shrink=0.7, label="f(x, y)")
    fig_sample.tight_layout()

    fig_uniform = plt.figure(figsize=(7, 5))
    ax_uniform = fig_uniform.add_subplot(111, projection="3d")
    scatter_uniform = ax_uniform.scatter(
        x_uniform,
        y_uniform,
        z_uniform,
        c=z_uniform,
        cmap=cmap,
        alpha=0.7,
        s=20,
    )
    ax_uniform.set_xlabel("x")
    ax_uniform.set_ylabel("y")
    ax_uniform.set_zlabel("f(x, y)")
    ax_uniform.set_xlim(*x_limits)
    ax_uniform.set_ylim(*y_limits)
    fig_uniform.colorbar(scatter_uniform, ax=ax_uniform, shrink=0.7, label="f(x, y)")
    fig_uniform.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
