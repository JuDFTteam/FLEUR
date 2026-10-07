---
name: fleur-bugfix
description: The FLEUR procedure for a bug found outside the current development scope - a pre-existing bug in the develop branch that is not part of what is being worked on right now. Opens a GitLab issue on iffgit.fz-juelich.de/fleur/fleur, records the fix plan and then the solution on that issue, commits the fix to develop with the issue referenced, and cherry-picks it onto the working branch. Use whenever a bug is found or reported in FLEUR and before starting to fix it, including when the bug turns up incidentally while working on something else.
---

# FLEUR bugfix procedure

Bugs that predate the current work belong to the whole FLEUR community, not to
one branch. This procedure makes sure such a fix is visible, traceable and
lands where everyone gets it: an issue that records the problem, the plan and
the solution, and a commit on `develop` that points back at that issue.

Work through the steps in order. Tell the user the issue URL as soon as it
exists, and the commit SHAs at the end.

## 1. Decide whether the procedure applies

It applies only if **all three** hold:

1. **Out of scope** — the bug is in code outside the direct scope of the current
   development effort. A bug in code being written right now is just part of
   that work: fix it and move on, no issue needed.
2. **Present in `develop`** — the faulty code is on the `develop` branch, not
   only on a feature branch.
3. **A real bug** — wrong results, a crash, a hang, a build break; not a known
   limitation, a work-in-progress gap or an unimplemented feature.

Verify gate 2 rather than guessing:

```bash
git fetch origin develop
git show origin/develop:src/fleur/<subdir>/<file>.F90 | sed -n '<first>,<last>p'   # is the faulty code there?
git blame -L <first>,<last> -- src/fleur/<subdir>/<file>.F90                       # who introduced it?
git branch -r --contains <blame-sha> | grep -x '  origin/develop'                  # is that commit on develop?
git log --oneline origin/develop..HEAD -- src/fleur/<subdir>/<file>.F90            # or did current work touch it?
```

If the blamed commit is already on `develop` and current work did not touch the
faulty lines, the procedure applies.

State in one sentence which gate decided it, then either continue or fix the bug
normally without the ceremony.

## 2. Get GitLab access sorted

Project: `fleur/fleur`, id **19**, on `https://iffgit.fz-juelich.de`. Reading
issues is public; creating issues and comments needs a token.

Use the bundled helper (stdlib Python only, no `glab` or `requests` needed):

```bash
SK=.claude/skills/fleur-bugfix/scripts/gitlab_issue.py
python3 $SK list --search "<keyword>"      # look for a duplicate first - always do this
python3 $SK show <iid>
python3 $SK create  --title "..." --body-file <file> --label Bug
python3 $SK comment <iid> --body-file <file>
```

It reads the token from `FLEUR_GITLAB_TOKEN`, `GITLAB_TOKEN`, `GITLAB_API_TOKEN`
or `CI_JOB_TOKEN` and exits with code **3** if none is set. If it exits 3
(and `glab` is not installed either), fall back to the manual route:

- write the issue body to a file in the scratchpad,
- show it to the user,
- point them at `https://iffgit.fz-juelich.de/fleur/fleur/-/issues/new` (the
  repo's `.gitlab/issue_templates/Bug.md` template is selectable there),
- ask for the issue number they got, and carry on from step 4 with it.

Opening an issue and commenting on it are visible to every FLEUR developer.
**Show the user the text and get their go-ahead before posting it** — each time,
unless they have said to stop asking. Never invent an issue number.

## 3. Open the issue — before writing the fix

Follow `.gitlab/issue_templates/Bug.md`, adapted for a developer-found bug.
Write it for someone who has never seen the code:

```markdown
# Summary

<One paragraph: what goes wrong, where. Name the subsystem and file:line.>

## Reproduction

<Test case, input file, or the code path that shows it. If it was found by
reading code rather than by running, say so and explain the trigger conditions.>

## This is a BUG because

<The expected behaviour, and why the present behaviour differs from it.>

## Implications

<Who is affected: which features (SOC, noco, DFPT, hybrid, LOs, film, ...),
serial vs. MPI vs. GPU, silently wrong numbers vs. crash. Since which commit or
release, if known. This is the part that tells other developers whether their
results are affected -- do not skip it.>

## The problem only occurs if

<Conditions that narrow it down, or "always".>

## Compute environment

<Compiler, machine, build flags -- if the bug is environment-dependent.>

## Ideas for fixes

<Initial hunch. The actual plan follows as a comment in step 4.>

/label ~Bug
```

Use `~"Critial Bug"` (the label really is spelled that way in this project)
instead of `~Bug` for wrong physics in a standard calculation, a data-corrupting
bug, or anything that makes published results untrustworthy.

Record the issue iid (`#NNN`) and give the user the URL.

## 4. Post the plan as a comment, as soon as it exists

Do this **before** implementing, as soon as the plan is settled — the issue
should show the reasoning, not just the outcome. Include:

- **Root cause** — the actual mechanism, not the symptom.
- **Files to change**, with one line each on what changes.
- **Approach**, and the alternatives rejected, with the reason.
- **Risk** — what else touches this code, which features could regress.
- **Test plan** — which existing tests cover it, which new or adapted test in
  `testing/tests/` will demonstrate the fix.
- **Cherry-pick needed?** — which branch also needs the fix.

If the plan changes while implementing, say so in the step 6 comment.

## 5. Implement and verify

Follow the conventions in `CLAUDE.md`: 3-space indent, `m_`-prefixed modules,
`implicit none` + `private`, `judft_error`/`judft_warn` instead of `stop`,
assumed-shape array dummies, no file I/O outside `io/`.

```bash
cd build && make -j                                   # or ./configure.sh -make
./run_tests.sh -k <substring>                         # the tests that cover the fix
./run_tests.sh -m <marker>                            # e.g. noco, soc, dfpt, forces
pre-commit run --files <changed files>                # copyright header, implicit none, no stop
fortitude check src/fleur/<subdir>/                   # static analysis
```

Add or adapt a regression test in `testing/tests/` when the bug is reachable
from a test input — a fix without a test invites the bug back. If no test is
practical, say why in the issue comment.

Keep the fix minimal and separate from current development work: only the files
the bug needs. Unrelated cleanups make the cherry-pick in step 8 fail.

## 6. Document the solution on the issue

After implementing, comment on the issue with:

- **Root cause**, confirmed or corrected against the step 4 plan.
- **What changed** — each file with a one-line description of the change.
- **Why this is the right fix**, and what was ruled out.
- **Test evidence** — the command run and its result, quoted, not paraphrased.
- **Remaining limitations** or follow-up work, if any.

Then post the commit SHA as a short follow-up comment once step 7 is done, so
the issue links to the code.

## 7. Commit on `develop`

The fix must originate on `develop` — that is the branch everyone builds from,
and the default branch of the project.

**If the current branch is `develop`:** stage only the bugfix files and commit.

**If the current branch is a feature branch** (the usual case), commit on
`develop` in a worktree, leaving the feature branch's working tree untouched:

```bash
WT=/tmp/fleur-bugfix-wt
git diff -- <bugfix files> > /tmp/fleur-bugfix.patch    # the fix, alone
git worktree add "$WT" develop
git -C "$WT" apply /tmp/fleur-bugfix.patch
git -C "$WT" add -- <bugfix files>
git -C "$WT" commit                                     # message below
git -C "$WT" log -1 --format=%H
```

Then drop the fix from the feature branch's working tree (`git checkout --
<bugfix files>`) so the cherry-pick in step 8 applies cleanly.

Commit message, matching the style already in this repo's history:

```
<Imperative subject, <= 72 chars, says what is fixed>

<Why the code was wrong and what the fix does. Wrap at 72 columns.
Mention affected features if results could change.>

Fixes #<iid>
```

`develop` is the default branch, so `Fixes #NNN` / `Closes #NNN` **auto-closes
the issue** when the commit lands there. Use it for a complete fix; use
`See #NNN` or `Partial fix for #NNN` when work remains. Append whatever
co-author attribution line the session is configured to add.

Pushing puts the fix in front of the whole project and starts a CI pipeline:

```bash
git -C "$WT" push origin develop     # only after the user says to
```

**Ask the user before pushing**, every time. The GitLab CI pipeline must pass;
if it breaks, fixing or reverting it is part of this job.

## 8. Cherry-pick onto the branch being worked on

Only if the current work actually needs the fix (it usually does — that is how
the bug was found):

```bash
git switch <feature-branch>
git cherry-pick <sha-from-step-7>
```

If the feature branch is meant to track `develop` anyway, `git merge develop`
is the better move — `CONTRIBUTING.md` asks for frequent merges. Either way,
rebuild and re-run the tests that the current work depends on.

Clean up the worktree when done:

```bash
git worktree remove /tmp/fleur-bugfix-wt
```

## Checklist

- [ ] Gate checked: out of scope, present in `develop`, a real bug
- [ ] Searched for an existing issue on the same bug
- [ ] Issue opened with summary, reproduction, **implications**, label
- [ ] Plan posted as a comment before implementing
- [ ] Fix implemented, built, tested; regression test added or its absence justified
- [ ] `pre-commit` clean on the changed files
- [ ] Solution documented on the issue
- [ ] Committed on `develop` with `Fixes #NNN`; pushed only after the user agreed
- [ ] Commit SHA posted to the issue
- [ ] Cherry-picked (or merged) onto the working branch; worktree removed
- [ ] User told: issue URL, commit SHA, cherry-pick SHA
