import os
from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache

output_dir = os.path.expanduser("~/allen_cache")
manifest_path = os.path.join(output_dir, "manifest.json")

cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)
sessions = cache.get_session_table()
print(f"Total sessions: {len(sessions)}")
print(sessions.head())