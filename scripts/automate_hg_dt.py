#!/usr/bin/env python3
import os
import glob
import json
import time
import requests
import subprocess

# Set this to your Jules API Key
API_KEY = os.environ.get("JULES_API_KEY")
if not API_KEY:
    print("Error: Please set the JULES_API_KEY environment variable.")
    print("Example: export JULES_API_KEY='your_api_key'")
    exit(1)

API_URL = "https://jules.googleapis.com/v1alpha"
SOURCE = "sources/github/AkeBoss-tech/HiC-TAD-Library"
WORK_ORDERS_DIR = "spec/work_orders"

def get_pr_branch(pr_url):
    """Uses GitHub CLI to fetch the branch name of a given PR URL."""
    try:
        res = subprocess.run(
            ["gh", "pr", "view", pr_url, "--json", "headRefName"], 
            capture_output=True, 
            text=True, 
            check=True
        )
        data = json.loads(res.stdout)
        return data.get("headRefName")
    except Exception as e:
        print(f"Error fetching PR info using gh CLI: {e}")
        if isinstance(e, subprocess.CalledProcessError):
            print(f"CLI Error Output: {e.stderr}")
        return None

def main():
    headers = {
         "X-Goog-Api-Key": API_KEY,
         "Content-Type": "application/json"
    }
    
    # 1. Glob work orders and sort them
    files = sorted(glob.glob(os.path.join(WORK_ORDERS_DIR, "*.md")))
    if not files:
        print(f"No work orders found in {WORK_ORDERS_DIR}.")
        return

    print(f"Found {len(files)} work orders to process.")

    # Starting branch: synthesis
    starting_branch = "synthesis"

    for file in files:
        print(f"\n{'='*50}")
        print(f"Processing {file}")
        print(f"Base branch: {starting_branch}")
        print(f"{'='*50}")
        
        with open(file, "r") as f:
            content = f.read()

        title = os.path.basename(file).replace('.md', '')
        prompt = (
            f"Please implement the following HG-DT work order. The spec details are below:\n\n{content}\n\n"
            "CRITICAL: After implementing the changes, you MUST run the appropriate test scripts or verification commands "
            "(e.g., 'make setup', 'make test', or equivalent) to ensure there are no regressions. "
            "You have access to all necessary environment variables in your workspace."
        )

        payload = {
            "prompt": prompt,
            "sourceContext": {
                "source": SOURCE,
                "githubRepoContext": {
                    "startingBranch": starting_branch
                }
            },
            "automationMode": "AUTO_CREATE_PR",
            "title": title
        }

        # 2. Start Jules session
        print("Creating Jules session...")
        try:
            resp = requests.post(f"{API_URL}/sessions", headers=headers, json=payload)
            resp.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"Failed to create session: {e}")
            if hasattr(e, 'response') and e.response is not None:
                print(e.response.text)
            break

        session = resp.json()
        session_id = session['id']
        session_name = session['name']
        print(f"Session Created successfully: {session_id}")

        # 3. Poll until complete and extract PR
        pr_url = None
        session_completed = False
        last_feedback_time = 0
        
        while not session_completed:
            time.sleep(15)
            print(f"Polling session {session_id} for completion...")
            
            try:
                poll_resp = requests.get(f"{API_URL}/{session_name}", headers=headers)
                
                if poll_resp.ok:
                    poll_session = poll_resp.json()
                    
                    if poll_session.get("state") == "AWAITING_USER_FEEDBACK":
                        current_time = time.time()
                        if current_time - last_feedback_time > 60:
                            print(f"Session {session_id} is awaiting feedback. Auto-confirming...")
                            msg_payload = {"prompt": "Yes, please proceed with the current plan. Make sure to run tests!"}
                            requests.post(f"{API_URL}/{session_name}:sendMessage", json=msg_payload, headers=headers)
                            last_feedback_time = current_time

                    if "outputs" in poll_session:
                        for output in poll_session["outputs"]:
                            if "pullRequest" in output:
                                pr_url = output["pullRequest"]["url"]
                                session_completed = True
                                break
                
                if not session_completed:
                    act_resp = requests.get(f"{API_URL}/{session_name}/activities?pageSize=10", headers=headers)
                    if act_resp.ok:
                        acts = act_resp.json().get("activities", [])
                        for activity in acts:
                            if "sessionCompleted" in activity:
                                session_completed = True
                                break
            except Exception as e:
                print(f"Error while polling: {e}")
                pass

        if pr_url:
            print(f"Session complete! PR created: {pr_url}")
            print("Extracting branch name from PR...")
            
            branch_name = get_pr_branch(pr_url)
            
            if branch_name:
                starting_branch = branch_name
                print(f"-> Next work order will build off branch: '{starting_branch}'")
            else:
                print("Could not determine PR branch. Aborting to avoid conflicts.")
                break
        else:
            print("Session completed, but no PR was produced. Stopping chain.")
            break

        print("Sleeping briefly before starting the next work order...")
        time.sleep(10)

    print("\nAll tasks processed! Check your repository for the generated PRs.")

if __name__ == '__main__':
    main()
