"""Standalone demo of the unified-input routing table.

Run with::

    streamlit run mock_routing_ui.py

Feeds hand-crafted RoutingDecision objects into the production UI function
so the layout and override flow can be evaluated without uploading anything.
Delete this file once the UI is wired into main.py.
"""

import streamlit as st

from utils.bloodgroup.router import CandidateScore, RoutingDecision
from utils.bloodgroup.router_ui import render_routing_table


st.set_page_config(layout="wide", page_title="Routing UI — mock")


def _mock_decisions() -> list[RoutingDecision]:
    """Four reads covering all three decision branches."""
    return [
        # Clear ABO match — auto-routed.
        RoutingDecision(
            read_id="patient_01_abo_exon7.fasta",
            sequence="N" * 720,
            decision="routed",
            best=CandidateScore("ABO", 98.7, "forward", 712),
            runner_up=CandidateScore("RHD", 71.2, "forward", 580),
            all_scores=[
                CandidateScore("ABO",  98.7, "forward", 712),
                CandidateScore("RHD",  71.2, "forward", 580),
                CandidateScore("RHCE", 70.8, "forward", 590),
                CandidateScore("KEL",  12.4, "reverse",  50),
            ],
            note="Routed to ABO at 98.7% (runner-up RHD 71.2%).",
        ),
        # Clear KEL match — auto-routed.
        RoutingDecision(
            read_id="patient_01_kel_exon6.ab1",
            sequence="N" * 540,
            decision="routed",
            best=CandidateScore("KEL", 96.4, "forward", 528),
            runner_up=CandidateScore("ABO", 68.1, "reverse", 200),
            all_scores=[
                CandidateScore("KEL",  96.4, "forward", 528),
                CandidateScore("ABO",  68.1, "reverse", 200),
                CandidateScore("RHD",  64.5, "forward", 180),
                CandidateScore("RHCE", 63.9, "forward", 175),
            ],
            note="Routed to KEL at 96.4% (runner-up ABO 68.1%).",
        ),
        # RHD/RHCE paralogy — ambiguous, needs review.
        RoutingDecision(
            read_id="patient_02_rh_amplicon.ab1",
            sequence="N" * 480,
            decision="ambiguous",
            best=CandidateScore("RHCE", 94.1, "forward", 460),
            runner_up=CandidateScore("RHD", 91.8, "forward", 458),
            all_scores=[
                CandidateScore("RHCE", 94.1, "forward", 460),
                CandidateScore("RHD",  91.8, "forward", 458),
                CandidateScore("ABO",  14.1, "reverse",  80),
                CandidateScore("KEL",  12.3, "forward",  45),
            ],
            note=("RHCE 94.1% vs RHD 91.8% (margin 2.3pp < 5.0pp) "
                  "— manual override recommended."),
        ),
        # Low-quality / off-target read — unrouteable.
        RoutingDecision(
            read_id="patient_03_unknown.fasta",
            sequence="N" * 300,
            decision="unknown",
            best=CandidateScore("RHD", 62.4, "forward", 180),
            runner_up=CandidateScore("ABO", 58.9, "reverse", 175),
            all_scores=[
                CandidateScore("RHD",  62.4, "forward", 180),
                CandidateScore("ABO",  58.9, "reverse", 175),
                CandidateScore("RHCE", 55.1, "forward", 170),
                CandidateScore("KEL",  41.8, "reverse",  90),
            ],
            note="Top identity 62.4% < 90.0% — no system matched confidently.",
        ),
    ]


st.title("🧪 Unified-input routing UI — mock")
st.caption(
    "Hand-crafted decisions, no real alignment. Adjust selections in the "
    "**Reads needing your review** section, then click confirm to see the "
    "final mapping the production code would feed into each analyzer."
)

decisions = _mock_decisions()
confirmed = render_routing_table(decisions)

if confirmed is not None:
    st.success("Routing confirmed!")
    st.markdown("#### Final mapping returned to caller")
    st.json({rid: (sys or "(skip)") for rid, sys in confirmed.items()})

    st.markdown("#### How the caller would dispatch this")
    by_system: dict[str, list[str]] = {}
    for rid, sys in confirmed.items():
        if sys is None:
            continue
        by_system.setdefault(sys, []).append(rid)
    for sys, rids in sorted(by_system.items()):
        st.write(f"- **{sys}** ← {len(rids)} read(s): {', '.join(rids)}")
    skipped = [rid for rid, sys in confirmed.items() if sys is None]
    if skipped:
        st.write(f"- _Skipped:_ {', '.join(skipped)}")
