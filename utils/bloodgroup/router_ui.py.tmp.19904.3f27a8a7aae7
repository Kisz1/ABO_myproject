"""Streamlit UI for the unified-input routing table.

Renders the routing decisions produced by `utils.bloodgroup.router.route_reads`,
provides per-row override controls for ambiguous / unknown reads, and returns
the user-confirmed routing mapping. Kept in a separate module from `router.py`
so the routing logic stays importable without pulling in Streamlit.
"""

from typing import Optional

import pandas as pd
import streamlit as st

from utils.bloodgroup.router import RoutingDecision


_STATUS_LABEL = {
    'routed':    '✓ Routed',
    'ambiguous': '⚠ Ambiguous',
    'unknown':   '✗ Unrouteable',
}


def _summary_row(d: RoutingDecision) -> dict:
    top = d.best
    runner = d.runner_up
    margin = (f"{top.identity - runner.identity:.1f}"
              if (top and runner) else '—')
    return {
        'Read': d.read_id,
        'Routed to': top.system if top else '—',
        'Identity (%)': f"{top.identity:.1f}" if top else '—',
        'Strand': top.strand if top else '—',
        'Runner-up': (f"{runner.system} {runner.identity:.1f}%"
                      if runner else '—'),
        'Margin (pp)': margin,
        'Status': _STATUS_LABEL[d.decision],
    }


def render_routing_table(
    decisions: list[RoutingDecision],
    confirm_label: str = "Confirm routing and run analysis",
) -> Optional[dict[str, Optional[str]]]:
    """Render the routing table plus override controls.

    Returns ``None`` while the user is still reviewing, and the final
    ``{read_id: system_or_None}`` mapping once they click the confirm
    button. ``None`` as a value means "skip this read" — the caller
    should drop it from the downstream analyzer inputs.
    """
    if not decisions:
        st.info("No reads to route.")
        return None

    counts = {'routed': 0, 'ambiguous': 0, 'unknown': 0}
    for d in decisions:
        counts[d.decision] += 1

    st.markdown("### 🧭 Routing Detection")
    c1, c2, c3 = st.columns(3)
    c1.metric("Auto-routed", counts['routed'])
    c2.metric("Ambiguous (needs review)", counts['ambiguous'])
    c3.metric("Unrouteable", counts['unknown'])

    st.dataframe(
        pd.DataFrame([_summary_row(d) for d in decisions]),
        hide_index=True,
        use_container_width=True,
    )

    # Per-row overrides for ambiguous / unknown reads. Auto-routed reads
    # are accepted as-is; user can still inspect them in the table above.
    overrides: dict[str, Optional[str]] = {}
    needs_review = [d for d in decisions if d.decision != 'routed']

    if needs_review:
        st.markdown("#### Reads needing your review")
        st.caption(
            "Pick the system to route each read to, or **Skip** to drop "
            "it from analysis. The suggested system is pre-selected."
        )
        for d in needs_review:
            st.markdown(f"**{d.read_id}** — {d.note}")

            score_lines = "\n".join(
                f"  - {s.system}: {s.identity:.1f}%  "
                f"({s.aligned_bases} aligned bases, {s.strand})"
                for s in d.all_scores
            ) or "  - (no scores)"
            st.code(score_lines, language=None)

            options = ['Skip'] + [s.system for s in d.all_scores]
            default_idx = (1 if d.best else 0)
            choice = st.selectbox(
                f"Route `{d.read_id}` to:",
                options,
                index=default_idx,
                key=f"_route_override_{d.read_id}",
            )
            overrides[d.read_id] = None if choice == 'Skip' else choice
            st.markdown("---")

    # Build the final mapping. Auto-routed reads use their winner;
    # reviewed reads use the user's selection (which may be None=Skip).
    final: dict[str, Optional[str]] = {}
    for d in decisions:
        if d.read_id in overrides:
            final[d.read_id] = overrides[d.read_id]
        elif d.decision == 'routed' and d.best:
            final[d.read_id] = d.best.system
        else:
            final[d.read_id] = None

    confirmed = st.button(confirm_label, type="primary")
    return final if confirmed else None
