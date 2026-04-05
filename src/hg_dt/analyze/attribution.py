from typing import Dict, Any

def generate_mechanistic_insight(mod_details: Dict[str, Any], delta_stats: Dict[str, Any]) -> str:
    """
    Generate a 'Mechanistic Insight' text string summarizing the multi-scale delta.

    Args:
        mod_details: Dictionary containing modification details like 'type', 'target', 'locus'.
        delta_stats: Dictionary containing computed deltas, e.g.,
                     {'accessibility_drop': 0.28, 'expression_drop': 0.35, 'loop_weakened': True}

    Returns:
        A human-readable string attributing the structural consequence.
    """
    mod_type = mod_details.get("type", "modification")
    target = mod_details.get("target", "element")

    insight_parts = [f"This {mod_type} affects {target}"]

    if delta_stats.get("loop_weakened"):
        insight_parts.append("weakens enhancer-promoter loop")
    elif delta_stats.get("loop_strengthened"):
        insight_parts.append("strengthens enhancer-promoter loop")

    acc_drop = delta_stats.get("accessibility_drop", 0)
    if acc_drop > 0:
        insight_parts.append(f"accessibility ↓{int(acc_drop * 100)}%")
    elif acc_drop < 0:
        insight_parts.append(f"accessibility ↑{int(-acc_drop * 100)}%")

    exp_drop = delta_stats.get("expression_drop", 0)
    if exp_drop > 0:
        insight_parts.append(f"expression ↓{int(exp_drop * 100)}%")
    elif exp_drop < 0:
        insight_parts.append(f"expression ↑{int(-exp_drop * 100)}%")

    # Example format: "This deletion removes ccRE EH38E1800647 -> weakens enhancer-promoter loop -> accessibility down 28% -> expression down 35%."
    return " → ".join(insight_parts) + "."
