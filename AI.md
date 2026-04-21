# AI Coding Assistants

This document provides guidance for AI tools and developers using AI assistance when contributing to `spyrmsd`.

AI tools helping with `spyrmsd` development should follow the standard development process.

## Licensing and Legal Requirements

All contributions must comply with the `spyrmsd`'s licensing requirements:

* All code must be compatible with the MIT license
* Use appropriate SPDX license identifiers


## Signed-off-by and Developer Certificate of Origin

AI agents MUST NOT add Signed-off-by tags. Only humans can legally certify the Developer Certificate of Origin (DCO). The human submitter is responsible for:

* Reviewing all AI-generated code
* Ensuring compliance with licensing requirements
* Adding their own Signed-off-by tag to certify the DCO
* Taking full responsibility for the contribution

## Attribution

When AI tools contribute to `spyrmsd` development, proper attribution helps track the evolving role of AI in the development process. Contributions should include an Assisted-by tag in the following format:

```
Assisted-by: AGENT_NAME:MODEL_VERSION [TOOL1] [TOOL2]
```

Where `AGENT_NAME` is the name of the AI tool or framework,
`MODEL_VERSION` is the specific model version used, and
`[TOOL1] [TOOL2]` are optional specialized analysis tools used.
Basic development tools (git, gcc, make, editors) should not be listed.
