"""Command line interface: ``sptpalm analyze <folder> --dt-ms 20 --pixel-size-um 0.102``."""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

from .params import AnalysisParams


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sptpalm", description="sptPALM analysis of TrackMate tracks")
    sub = parser.add_subparsers(dest="command", required=True)

    a = sub.add_parser("analyze", help="load the TrackMate spot tables of a folder and analyse them")
    a.add_argument("folder", help="folder with the TrackMate spots tables (one file per movie)")
    a.add_argument("--dt-ms", type=float, required=True, help="acquisition time (ms)")
    a.add_argument("--pixel-size-um", type=float, required=True, help="pixel size (um)")
    a.add_argument("--pattern", default="*.csv", help="file pattern of the tables (default *.csv)")
    a.add_argument("--out", default=None, help="output folder (default: <folder>/sptpalm_results)")
    a.add_argument("--max-blinks", type=int, default=3, help="maximum gap (frames) inside a track")
    a.add_argument("--max-step-um", type=float, default=None,
                   help="maximum step length (um); longer steps are treated as mis-connections")
    a.add_argument("--min-points", type=int, default=8, help="minimum number of detections per track")
    a.add_argument("--min-duration-frames", type=int, default=0)
    a.add_argument("--min-fraction", type=float, default=0.0,
                   help="minimum fraction of frames with a detection (0.75 recommended if used)")
    a.add_argument("--min-points-msd", type=int, default=3)
    a.add_argument("--msd-points", type=int, default=4, help="MSD points used to calculate D")
    a.add_argument("--max-display-time", type=float, default=0.5, help="MSD plotted up to this time (s)")
    a.add_argument("--method", choices=["average", "fit"], default="average")
    a.add_argument("--gaussians", type=int, choices=[1, 2], default=2)
    a.add_argument("--split", type=float, default=None,
                   help="approximate log10(D) between the two populations (default: automatic)")
    a.add_argument("--legacy", action="store_true", help="reproduce the MATLAB behaviour")
    a.add_argument("--no-plots", action="store_true", help="do not save the figures")
    a.add_argument("--trajectories", action="store_true", help="also plot the trajectories")
    return parser


def _analyze(args) -> int:
    from .io.export import write_results
    from .io.trackmate import load_trackmate_folder
    from .pipeline import run_analysis

    params = AnalysisParams(
        acquisition_time_ms=args.dt_ms, pixel_size_um=args.pixel_size_um,
        max_blinks=args.max_blinks, max_step_length_um=args.max_step_um,
        min_points=args.min_points, min_duration_frames=args.min_duration_frames,
        min_fraction=args.min_fraction, min_points_msd=args.min_points_msd,
        msd_fit_points=args.msd_points, max_display_time_s=args.max_display_time,
        diffusion_method=args.method, n_gaussians=args.gaussians, split_threshold=args.split,
        legacy=args.legacy)

    print(f"Loading {args.folder} ...")
    tracks, info = load_trackmate_folder(args.folder, params.pixel_size_um, args.pattern)
    if tracks.n_tracks == 0:
        print("No TrackMate track found.", file=sys.stderr)
        return 1
    print(f"{len(info.files)} file(s), {tracks.n_tracks} tracks "
          f"({info.n_rejected_total} rejected because of a zero-length step)")

    results = run_analysis(tracks, params)
    out = args.out or os.path.join(args.folder, "sptpalm_results")
    write_results(results, out)
    if not args.no_plots:
        from .plotting import save_figures
        save_figures(results, out, loaded_tracks=tracks, trajectories=args.trajectories)

    for name, r in results.items():
        print(f"[{name}] loaded {r.n_loaded} -> filtered {r.n_filtered} -> MSD {r.n_msd} -> D {r.n_dapp}")
        if r.fit is not None:
            for k, c in enumerate(r.fit.centers, start=1):
                print(f"    population {k}: log10(D) = {c:.3f}  (D = {10 ** c:.4f} um^2/s)")
            if r.fit.fast_fraction is not None:
                print(f"    mobile fraction: {100 * r.fit.fast_fraction:.1f}%")
        else:
            print(f"    no fit: {r.fit_error}")
    print(f"Results written in {out}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "analyze":
        return _analyze(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
