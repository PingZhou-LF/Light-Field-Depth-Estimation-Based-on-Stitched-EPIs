; Listing generated by Microsoft (R) Optimizing Compiler Version 19.28.29335.0 

include listing.inc

INCLUDELIB MSVCRT
INCLUDELIB OLDNAMES

msvcjmc	SEGMENT
__0C77F790_corecrt_stdio_config@h DB 01H
__ABC2E6F9_corecrt_wstdio@h DB 01H
__59E0AD71_stdio@h DB 01H
__561A04E6_stdlib@h DB 01H
__48296EF4_tmwtypes@h DB 01H
__E637FCD0_cpp_mexapi_version@cpp DB 01H
msvcjmc	ENDS
PUBLIC	mexfilerequiredapiversion
PUBLIC	__JustMyCode_Default
EXTRN	_RTC_InitBase:PROC
EXTRN	_RTC_Shutdown:PROC
EXTRN	__CheckForDebuggerJustMyCode:PROC
pdata	SEGMENT
$pdata$mexfilerequiredapiversion DD imagerel $LN3
	DD	imagerel $LN3+75
	DD	imagerel $unwind$mexfilerequiredapiversion
pdata	ENDS
;	COMDAT rtc$TMZ
rtc$TMZ	SEGMENT
_RTC_Shutdown.rtc$TMZ DQ FLAT:_RTC_Shutdown
rtc$TMZ	ENDS
;	COMDAT rtc$IMZ
rtc$IMZ	SEGMENT
_RTC_InitBase.rtc$IMZ DQ FLAT:_RTC_InitBase
rtc$IMZ	ENDS
xdata	SEGMENT
$unwind$mexfilerequiredapiversion DD 022301H
	DD	0700b320fH
xdata	ENDS
; Function compile flags: /Odt
;	COMDAT __JustMyCode_Default
_TEXT	SEGMENT
__JustMyCode_Default PROC				; COMDAT
	ret	0
__JustMyCode_Default ENDP
_TEXT	ENDS
; Function compile flags: /Odtp /RTCsu
_TEXT	SEGMENT
built_by_rel$ = 48
target_api_ver$ = 56
mexfilerequiredapiversion PROC
; File E:\matlab2020\matlab\extern\version\cpp_mexapi_version.cpp
; Line 9
$LN3:
	mov	QWORD PTR [rsp+16], rdx
	mov	QWORD PTR [rsp+8], rcx
	push	rdi
	sub	rsp, 32					; 00000020H
	mov	rdi, rsp
	mov	ecx, 8
	mov	eax, -858993460				; ccccccccH
	rep stosd
	mov	rcx, QWORD PTR [rsp+48]
	lea	rcx, OFFSET FLAT:__E637FCD0_cpp_mexapi_version@cpp
	call	__CheckForDebuggerJustMyCode
; Line 10
	mov	rax, QWORD PTR built_by_rel$[rsp]
	mov	DWORD PTR [rax], 131594			; 0002020aH
; Line 11
	mov	rax, QWORD PTR target_api_ver$[rsp]
	mov	DWORD PTR [rax], 120586240		; 07300000H
; Line 12
	add	rsp, 32					; 00000020H
	pop	rdi
	ret	0
mexfilerequiredapiversion ENDP
_TEXT	ENDS
END