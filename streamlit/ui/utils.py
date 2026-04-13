from pathlib import Path

def detect_group(p: Path, root: Path) -> str:
    """
    Detects the experimental group (CTL or Hypoxia) of a file based on its path.
    """
    try:
        rel = str(p.relative_to(root)).lower()
    except Exception:
        rel = str(p).lower()
    
    if "/hip/" in rel:
        return "Hypoxia"
    if "/ctl/" in rel:
        return "CTL"
    
    # Default to CTL if no specific group folder is found
    return "CTL"
